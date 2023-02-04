#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <sstream>

#include <cstdlib>
#include <cstring>
#include <htslib/sam.h>

#include "Tokenizer.h"

using namespace std;

int ReadLength = 0;

class Node
{
	public:

	string nodeid;
	
	// Chromosome.
	string chromosome;

	// Starting and ending coordinate.
	int start, end;

	// Forward or reverse strand.
	bool strand;

	// Ids.
	set<string> geneid, transcriptid;
	
	Node(string nodeid, string chromosome, int start, int end, bool strand) :
		nodeid(nodeid), chromosome(chromosome), start(start), end(end), strand(strand)
	{ }

	Node(const Node& N) :
		nodeid(N.nodeid), chromosome(N.chromosome), start(N.start), end(N.end), strand(N.strand)
	{ }	

	Node() { }

	// The amount of bases in the subexon.
	int size()
	{
		return end - start + 1;
	}

    bool contains(int x)
    {
        return start >= x && x <= end;
    }
	
	// Some operators for easier comparison.
	bool operator<(const Node& rhs) const
	{
		return start < rhs.start;
	}
	
	bool operator<=(const Node& rhs) const
	{
		return start <= rhs.start;
	}
	
	bool operator==(const Node& rhs) const
	{
		return start == rhs.start;
	}
};

// Data container for the GTF transcript
class Transcript
{
	public:

	// Subexons belonging to a transcript.
	set<Node*> nodes;

	// Transcriptomic to genomic coordinates.
	// x - starting position within the transcript.
	// rev - is transcript on reverse strand?
	// trim - is transcript trimmed?
	int t2g(int x, bool rev, bool trim)
	{
		// Sorting nodes by their starting coordinates.
		vector<Node*> v(nodes.begin(), nodes.end());
		sort(v.begin(), v.end(), [](const Node* l, const Node* r) 
		{ 
			return *l < *r; 
		});
		
        // Trim left/right.
        trim = trim && v.size() > 1;
        bool tl = trim && v.front()->size() >= ReadLength; 
        bool tr = trim && v.back()->size() >= ReadLength;

		// Offseting start for reverse strand.
        if (rev) 
        {
            int c = 0;
            for (int i = 0; i < (int)v.size(); ++i)
			{
                c += v[i]->size();
			}
            x = c - (tl + tr + 1) * (ReadLength - 1) - (x - 1);		
        }

		// Applying trimming criteria if applicable.
		if (tl) 
		{
			x += v.front()->size() - ReadLength + 1;
		}
	
		// Iterates through subexons in a transcript, finding to which one
		// does starting position x belong to.
		for (int i = 0; i < (int)v.size(); ++i)
		{
			if (x <= v[i]->size()) 
			{
				return max(-1, v[i]->start + x - 1);
			}
            else
			{
			    x -= v[i]->size();
			}
		}

		// x not found, default to -1.
		return -1;
	}

    int g2t(int x)
    {
        vector<Node*> v(nodes.begin(), nodes.end());

		sort(v.begin(), v.end(), [](const Node* l, const Node* r) 
		{ 
			return *l < *r; 
		});

        int c = 0;

        for (int i = 0; i < (int)v.size(); ++i)
		{
            if (v[i]->start <= x && x <= v[i]->end)
			{
                return c + x - v[i]->start + 1;
			}
            else
			{
                c += v[i]->size();
			}
		}

        return -1;
    }
    
};

// GTF parser/container.
class GTF
{
	public:

	GTF(string filename)
	{
		in.open(filename);
		while(getLine()) parseLine();
		line.clear();		
		in.close();			
	}

    bool nodeExists(string nodeid)
    {
        return N.find(nodeid) != N.end();
    }

	// Returns node for the specified id. Assumes that the node exists.
	Node& getNode(string nodeid)
	{
		return N[nodeid];
	}

	// Returns transcript for the specified id. Assumes that the transcript exists.
	Transcript& getTranscript(string transcriptid)
	{
		return T[transcriptid];
	}

	private:

	bool getLine()
	{
		for (string l; getline(in, l); )
		{
			l = l.substr(0, l.find("#"));			

			if (!l.size()) 
			{
				continue;
			}

			line.clear();
			line.str(l);
			return true;
		}
		return false;
	}

	bool nextToken()
	{
		if (line >> token)
		{
			if (token.back() == ';')
			{
				token = token.substr(0, token.size() - 1);
			}

			if (token.back() == '\"')
			{
				token = token.substr(1, token.size() - 2);
			}

			return true;
		}
		return false;
	}

	void parseLine()
	{
		int start = -1; 
		int end = -1;
		bool strand = 0;
		string chromosome = "";
		string nodeid = "";
		string geneid = ""; 
		string transcriptid = "";
		
		for(int i = 0; nextToken(); ++i)
		{
			if (i == 0)
			{
				chromosome = token;
			}
			else if (i == 2 && token != "subexon")
			{
				return;	
			}
			else if (i == 3)
			{
				start = stoi(token);
			}
			else if (i == 4)
			{
				end = stoi(token);
			}
			else if (i == 6)
			{
				strand = token == "+";
			}
			else if (token == "NodeId")
			{
				nextToken(); 
				nodeid = token;
			}
			else if (token == "gene_id")
			{
				nextToken();
				geneid = token;
			}
			else if (token == "transcript_id")
			{
				nextToken();
				transcriptid = token;
			}
		}

		Node& R = N[nodeid];
		R.nodeid = nodeid;
		R.chromosome = chromosome;
		R.start = start;
		R.end = end;
		R.strand = strand;
		R.geneid.insert(geneid);
		R.transcriptid.insert(transcriptid);
		T[transcriptid].nodes.insert(&R);
	}

	map<string, Node> N;
	map<string, Transcript> T;
	
	string token;
	istringstream line;
	ifstream in;
};

class BAM
{
	public:

	BAM(string filename)
	{
		IN = hts_open(filename.c_str(), "r");
		HD = sam_hdr_read(IN);
		A  = bam_init1();
	}

	bool next()
	{
		return sam_read1(IN, HD, A) > 0;
	}

	~BAM()
	{
		bam_destroy1(A);
		sam_close(IN);
	}

	samFile*   IN;
	bam_hdr_t* HD;
	bam1_t*    A;
};

vector<string> SplitString(string s, char c)
{
  	vector<string> V({""});

	for(unsigned int i = 0; i < s.size(); i++)
	{
		if (s[i] != c)
		{
			V.back() += s[i];
		}
		else if (i < s.size() - 1)
		{
			V.push_back("");
		}
	}
		
	return V;
}

Transcript AS_Frag2Trans(GTF& G, vector<string>& V)
{
	Transcript T;

	for (unsigned int i = 0; i < V.size(); ++i)
	{
		T.nodes.insert(&G.getNode(V[i]));
	}

	return T;
}

string AS_Frag2Cigar(Transcript& T, int s)
{
	vector<Node*> v(T.nodes.begin(), T.nodes.end());
	s = T.t2g(s, false, true);
	string cig = "";

	sort(v.begin(), v.end(), [](const Node* l, const Node* r) 
	{ 
		return *l < *r; 
	});

	for (unsigned int i = 0, k = 0, c = 0; i < v.size(); ++i)
	{
		if (v[i]->end < s)
		{
			continue;
		}

		int t = min(ReadLength - (int)k, v[i]->end - max(v[i]->start, s) + 1);
		int l = i < v.size() - 1 ? v[i + 1]->start - v[i]->end - 1 : 0; 
		k += t;
		c += t;
		
		bool f = i == v.size() - 1 || (int)k >= ReadLength;

		if (l || f)
		{
			cig += to_string(c) + "M" + (f ? "" : to_string(l) + "N");
			c = 0;		
		}	

		if ((int)k >= ReadLength)
		{
			return (cig.size() && cig.back() == 'N' ? "" : cig);
		}
	}

	return "";
}

string conv(unsigned char* A, int n)
{
	string s(n, ' ');

	for (int i = 0; i < n; ++i)
	{
		switch (bam_seqi(A, i))
		{
			case 1:  s[i] = 'A'; break;
			case 2:  s[i] = 'C'; break;
			case 4:  s[i] = 'G'; break;
			case 8:  s[i] = 'T'; break;
			default: s[i] = 'N';
		}
	}

	return s;
}

void P2G_All(string bam, string gtf, string fa)
{
	cerr << "Loading GTF file... " << flush;
	GTF G(gtf); BAM B(bam);
	cerr << "done. " << endl;
	
	cerr << "Loading FA file... " << flush;
	Tokenizer F(fa);
	vector<string> cn;
	vector<int> ci;
	map<string, int> tid;

	while (F.nextLine())
	{
		string s = F.getToken();

		if (!s.size())
		{
			continue;
		}

		if (s[0] == '>')
		{
			s = s.substr(1, s.size());
			cn.push_back(s);
			ci.push_back(0);
		}
		else
		{
			ci.back() += s.size();
		}
	}
	cerr << "done. " << endl;

	cerr << "Processing BAM file... " << flush;
  	std::string text = "@HD\tVN:1.0\n@PG\tID:kallisto\tPN:kallisto\tVN:1.0\n";
	int num_chr = cn.size();
    for (int i = 0; i < num_chr; i++)
	{
    	text += "@SQ\tSN:" + cn[i] + "\tLN:" + to_string(ci[i]) + "\n";
	}
	cout << text << flush;

	int wrongLen = 0;
	int wrongCig = 0;

	while (B.next())
	{
        if (B.A->core.l_qseq != ReadLength) 
		{ 
			++wrongLen; 
			continue;
		}

		if (B.A->core.flag == 4)
		{
			continue;
		}
	
		string frag = B.HD->target_name[B.A->core.tid];
		vector<string> V = SplitString(frag, '|');
		string fa = conv(bam_get_seq(B.A), ReadLength);
		Transcript T = AS_Frag2Trans(G, V);
		string cig = AS_Frag2Cigar(T, B.A->core.pos + 1);
		
		if (cig == "") 
		{
			++wrongCig; 
			continue;
		}

		cout << bam_get_qname(B.A) 			          << '\t'
             << B.A->core.flag      			      << '\t'
		     << (*T.nodes.begin())->chromosome        << '\t'
			 << T.t2g(B.A->core.pos + 1, false, true) << '\t'
			 << "255"			    			      << '\t'
			 << cig                                   << '\t'
			 << "*\t0\t0"                             << '\t'
			 << fa                                    << '\t'
			 << string(ReadLength, 'I')               << '\t'
			 << bam_aux2i(bam_aux_get(B.A, string("NH").c_str())) << endl;
	}
	cerr << " done." << endl;
	
	cerr << "wrong cig: " << wrongCig << '\n'
         << "incorrect len: " << wrongLen << endl;
}

int main(int argc, char** argv)
{

	if (argc != 9 && argc != 5)
	{
		cerr << "Converts a fortuna pseudobam into a genomesam.\n"
             << "Usage: ./run <bam> <gtf> <fa> <read length>" 
			 << endl;
	}
	else if (argc == 5)
	{
		ReadLength = stoi(argv[4]);
		P2G_All(argv[1], argv[2], argv[3]);
	}

	return 0;
}
