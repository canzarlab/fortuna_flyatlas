#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "Tokenizer.h"

using namespace std;

map<string, string> CellrangerMap(string cellrangerSamPath)
{
	Tokenizer T(cellrangerSamPath);
	map<string, string> M;

    while (T.nextLine())
    {
        if (T.getLine()[0] == '@')
        {
        	continue;
        }

		string read = T.getToken();

        for (int i = 0; i < 10; ++i)
        {
            T.getToken();
        }
                
        string token = T.getToken();

		string barcode = "";
		int flags = 0;
		bool readFlag = false;
		bool readBarcode = false;
        while (token != "")
        {
            string prefix = token.substr(0, 5);

            if (prefix == "xf:i:")
            {
                flags = stoi(token.substr(5, token.size() - 5));
				readFlag = true;
            }
			else if (prefix == "CB:Z:")
			{
				unsigned int index = token.find('-');
				if (index == string::npos)
				{
					barcode = token.substr(5, token.size() - 5);
				}
				else
				{
					barcode = token.substr(5, index - 5);
				}
				
				readBarcode = true;
			}

			if (readFlag && readBarcode)
			{
				break;
			}

            token = T.getToken();
        }

		if (readFlag && readBarcode && (flags & 8))
		{
			M[read] = barcode;
		}

    }
	
	return M;
}

map<tuple<string, int, int>, vector<string>> FortunaAS(string fortunaASPath)
{
	Tokenizer T(fortunaASPath);
	map<tuple<string, int, int>, vector<string>> M;

    while (T.nextLine())
    {
		string chromosome = T.getToken();
		string start = T.getToken();
		string end = T.getToken();
                
		// skip 1 value
        T.getToken();

		string type = T.getToken();

		tuple<string, int, int> key(chromosome, stoi(start), stoi(end));

		auto it = M.find(key);
		if (it == M.end())
		{
			vector<string> vec;
			vec.push_back(type);
			M[key] = vec;
		}
		else
		{
			it->second.push_back(type);
		}

    }
	
	return M;
}

map<tuple<string, int, int>, vector<string>> Count(map<string, string>& CMap, map<tuple<string, int, int>, vector<string>>& fortunaASMap, string fortunaSamPath)
{
	Tokenizer T(fortunaSamPath);
	map<tuple<string, int, int>, vector<string>> M;

    while (T.nextLine())
    {
        if (T.getLine()[0] == '@')
        {
        	continue;
        }

		string read = T.getToken();
		auto CMapIt = CMap.find(read);
		if (CMapIt == CMap.end())
		{
			continue;
		}

		// skip 1 val		
		T.getToken();

		string chromosome = T.getToken();
		string start = T.getToken();
		int startCoord = stoi(start);
		
		// skip 1 val		
		T.getToken();
		
		string cigar = T.getToken();

		// parse cigar string
		int prevNum = 0;
		unsigned int startIndex = 0;
		while (startIndex < cigar.length())
		{
			unsigned int endIndex = cigar.find_first_of("MNS", startIndex);
			char c = cigar[endIndex];

			if (c == 'S')
			{
				continue;
			}

			int num = stoi(cigar.substr(startIndex, endIndex - startIndex));

			if (c == 'N')
			{
				int intronStart = startCoord + prevNum - 1;
				int intronEnd = intronStart + num + 1;

				tuple<string, int, int> key(chromosome, intronStart, intronEnd);
				if (fortunaASMap.find(key) != fortunaASMap.end())
				{
					if (M.find(key) == M.end())
					{
						M[key] = vector<string>({ CMapIt->second });
					}
					else
					{
						M[key].push_back(CMapIt->second);
					}
				}

				startCoord += num + prevNum;
			}
			
			prevNum = num;
			startIndex = endIndex + 1;
		}
    }

	return M;
}

void GenerateAS(string fortunaSamPath, string cellrangerSamPath, string fortunaASPath)
{
    map<string, string> CMap = CellrangerMap(cellrangerSamPath);
	map<tuple<string, int, int>, vector<string>> fortunaASMap = FortunaAS(fortunaASPath);
	map<tuple<string, int, int>, vector<string>> countMap = Count(CMap, fortunaASMap, fortunaSamPath);

	for (auto& it : countMap)
	{
		vector<string>& events = fortunaASMap[it.first];

		cout << get<0>(it.first) << '\t' << get<1>(it.first) << '\t' << get<2>(it.first) << '\t';
        
		for (unsigned int i = 0; i < events.size(); ++i)
		{
			cout << events[i] << (i == events.size() - 1 ? '\t' : ',');
		}

		cout << it.second.size() << '\t';

		for (unsigned int i = 0; i < it.second.size(); ++i)
		{
			cout << it.second[i] << (i == it.second.size() - 1 ? '\n' : ',');
		}
	}
	cout << flush;
}

int main(int argc, char** argv)
{
    if (argc == 4)
    {
        GenerateAS(argv[1], argv[2], argv[3]);
    }
    else
    {
        cerr << "Adjusts fortuna’s catalog of novel AS events according to the information obtained from fortuna's pseudoalignment and cellranger's alignment.\n"
             << "Usage: ./run <fortuna-GENOME-SAM> <cellranger-SAM> <fortuna-ALT-SPLICING-TSV>" 
             << endl;
    }

    return 0;
}
