#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

class Tokenizer
{
	public:

	Tokenizer(string filename)
	{
		in.open(filename);		
	}

	string getToken()
	{
        return (ss >> token) ? token : "";
	}

	string getLine()
	{
		return line;
	}

    bool nextLine()
    {
        while(getline(in, line))
		{			
			if (!line.size()) 
			{
				continue;
			}

			ss.clear();
			ss.str(line);

			return true;
		}

		return false;
    }

    private:

	string token, line;
	istringstream ss;
	ifstream in;
};
