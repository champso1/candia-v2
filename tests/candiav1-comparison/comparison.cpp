#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <exception>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
using namespace std;


#define ANSI_COLOR_RED         "\x1b[31m"
#define ANSI_COLOR_GREEN       "\x1b[32m"
#define ANSI_COLOR_YELLOW      "\x1b[33m"
#define ANSI_COLOR_BLUE        "\x1b[34m"
#define ANSI_COLOR_MAGENTA     "\x1b[35m"
#define ANSI_COLOR_CYAN        "\x1b[36m"

#define ANSI_COLOR_RESET       "\x1b[0m"



using uint = unsigned int;

static string const& output_prefix = "output";
static string const& v1_output_prefix = "output-candiav1";
static string const& v2_output_prefix = "output-candiav2";

using FileContents = vector<vector<double>>;
vector<string> SplitString(string const& str, string const& delim);
FileContents GrabFile(string const& filepath);
FileContents GrabFile2(string const& filepath, const uint line_len);
void CompareFile(FileContents v1, FileContents v2, bool print);

int main()
{
	cerr << "testing LO...\n";
	FileContents LO_v1   = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-LO-splitfuncs.dat", 16);
	FileContents LO_v2   = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-LO-splitfuncs.dat", 16);
	CompareFile(LO_v1, LO_v2, false);
	cerr << "LO okay\n";

	cerr << "testing NLO...\n";
	FileContents NLO_v1  = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-NLO-splitfuncs.dat", 19);
	FileContents NLO_v2  = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-NLO-splitfuncs.dat", 19);
	CompareFile(NLO_v1, NLO_v2, false);
	cerr << "NLO okay\n";

	cerr << "testing NNLO...\n";
	FileContents NNLO_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-NNLO-splitfuncs.dat", 22);
	FileContents NNLO_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-NNLO-splitfuncs.dat", 22);
	CompareFile(NNLO_v1, NNLO_v2, true);
	cerr << "NNLO okay\n";

	cerr << "testing alphas_betas...\n";
	FileContents alphas_betas_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-alphas_betas.dat", 6);
	FileContents alphas_betas_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-alphas_betas.dat", 6);
	CompareFile(alphas_betas_v1, alphas_betas_v2, false);
	cerr << "Alphas and betas okay\n";

	cerr << "testing init_dists...\n";
	FileContents init_dists_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-init_dists.dat", 7);
	FileContents init_dists_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-init_dists.dat", 7);
	CompareFile(init_dists_v1, init_dists_v2, false);
	cerr << "initial distributions okay\n";
	
	cerr << "testing gauleg...\n";
	FileContents gauleg_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-gauleg.dat", 2);
	FileContents gauleg_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-gauleg.dat", 2);
	CompareFile(gauleg_v1, gauleg_v2, false);
	cerr << "gauleg stuff okay\n";

	cerr << "testing interp...\n";
	FileContents interp_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-interp.dat", 1);
	FileContents interp_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-interp.dat", 1);
	CompareFile(interp_v1, interp_v2, false);
	cerr << "interp stuff okay\n";

	cerr << "testing conv...\n";
	FileContents conv_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "out-conv.dat", 1);
	FileContents conv_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "out-conv.dat", 1);
	CompareFile(conv_v1, conv_v2, false);
	cerr << "conv stuff okay\n";

	cerr << "testing actual coefficient values from main evolution...\n";
	FileContents A_1_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "1.dat", 801);
	FileContents A_1_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "1.dat", 801);
	FileContents A_2_v1 = GrabFile2(output_prefix + '/' + v1_output_prefix + '/' + "2.dat", 801);
	FileContents A_2_v2 = GrabFile2(output_prefix + '/' + v2_output_prefix + '/' + "2.dat", 801);
	CompareFile(A_1_v1, A_1_v2, false);
	cerr << "coefficients from first iteration okay\n";
	CompareFile(A_2_v1, A_2_v2, false);
	cerr << "coefficients from second iteration okay\n";
	cerr << "all coefficients okay\n";
	

	return 0;
}


vector<string> SplitString(string const& _str, string const& _delim)
{
	const char* str = _str.c_str();
	const char* delim = _delim.c_str();
	vector<string> res;

	// copy the string data into a mutable buffer
	// so that we can tokenize it
	char buf[256];
	memset(buf, 0, 256);
	strncpy(buf, str, (sizeof buf)-1);

	char* ptr;
	ptr = strtok(buf, delim);
	while (ptr != nullptr) {
		res.emplace_back(ptr);
		ptr = strtok(nullptr, delim);
	}

	return res;
}


FileContents GrabFile(string const& filepath)
{
	ifstream f(filepath);
	if (!f)
	{
		cerr << "Could not open file " << filepath << '\n';
		exit(1);
	}

	FileContents res;

	string temp;
	while (!f.eof()) {
		getline(f, temp);
		vector<string> tokens = SplitString(temp, "\t");
		vector<double> vals;
		for (string const& tok : tokens)
		{
			try
			{
				vals.push_back(stod(tok));
			}
			catch (exception& e)
			{
			    continue;
			}
		}
		res.emplace_back(vals);
	}
	f.close();

	return res;
}

FileContents GrabFile2(string const& filepath, const uint line_len)
{
	ifstream f(filepath);
	if (!f)
	{
		cerr << "Error opening file '" << filepath << "'\n";
		exit(1);
	}

	FileContents contents;
	string temp;
	
	while (getline(f, temp))
	{
		if (temp.size() <= 1)
			continue;
		
		istringstream stream(temp);
		std::vector<double> data;
		string token;

		while (getline(stream, token, '\t'))
		{
			try
			{
				data.push_back(stold(token));
			}
			catch (invalid_argument const&)
			{
				continue;
			}
		}

		if (data.size() != line_len)
		{
			cerr << "Number of entries in file '" << filepath << "' is " << data.size()
				 << ", not the expected " << line_len << ".\n";
			exit(1);
		}

		contents.emplace_back(data);
	}
    
	f.close();
	return contents;
}


void CompareFile(FileContents v1, FileContents v2, bool print)
{
	(void)print;
	
	if (v1.size() != v2.size())
	{
		cerr << "Number of lines differ! " << v1.size() << "!=" << v2.size() << '\n';
		exit(1);
	}
	if (v1[0].size() != v2[0].size())
	{
		cerr << "Number of entries differ! " << v1[0].size() << "!=" << v2[0].size() << '\n';
		exit(1);
	}

	for (uint i=0; i<v1.size(); i++)
	{
		for (uint j=0; j<v1[i].size(); j++)
		{
			double avg = (v1[i][j] + v2[i][j])/2.0;
			double dif = abs(v1[i][j] - v2[i][j]);
			double dif_percentage = dif*100.0/avg;
			
			if (abs(dif_percentage) > 1e-2)
			{
				char x;
				cerr << "Values at position [" << i << ',' << j << "] = ["
					 << v1[i][j] << ',' << v2[i][j] << "] differ by more than .01%\n"
					 << "Continue? (Y/n) ";
				cin >> x;

				if (x == 'Y' || x == 'y')
					continue;
				else
					exit(1);
			}
		}
	}
}

