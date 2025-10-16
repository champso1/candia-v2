#include <cstdlib>
#include <ios>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <print>
using namespace std;

#include "Candia-v2/Candia.hpp"
using namespace Candia2;

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
	    println("[ERROR] Must provide a file to read from.");
		exit(EXIT_FAILURE);
	}

	string filename{argv[1]};

	filesystem::path filepath{filename};
	if (!filesystem::exists(filepath))
	{
	    println("[ERROR] Failed to find file \"{}\"", filepath.string());
		exit(EXIT_FAILURE);
	}
	
	ifstream file_stream{filepath};
	file_stream.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore the comment line

	vector<double> xtab{};
	vector<int> ntab{};
	
	double temp1{};
	int temp2{};
	string line{};

	// read in xtab array
	getline(file_stream, line);
	istringstream iss{line};
	while (iss >> temp1)
		xtab.push_back(temp1);

	// read in ntab array
	getline(file_stream, line);
	iss = istringstream{line};
	while (iss >> temp2)
		ntab.push_back(temp2);

	// read in rest of data points
	vector<double> X{};
	MultiDimVector<double, 2>::type F{};
	F.resize(13);
	while (getline(file_stream, line))
	{
		iss = istringstream{line};
		iss >> temp1;
		X.push_back(temp1);
		for (int i=0; i<F.size(); ++i)
		{
			iss >> temp1;
			F.at(i).push_back(temp1);
		}
	}

	// must create required dists from the regular ones
	MultiDimVector<double, 2>::type dists(8, vector<double>(11, 0.0));
	double size = ntab.size();
	for (uint ik=0; ik<size-1; ++ik) // -1 because we don't want x=1.0x
	{
		int k = ntab.at(ik);
		
		dists.at(0).at(ik) = F[1][k] - F[1+6][k];
		dists.at(1).at(ik) = F[2][k] - F[2+6][k];
		dists.at(2).at(ik) = F[2+6][k] - F[1+6][k];
		dists.at(3).at(ik) = 2.0*(F[2+6][k] + F[1+6][k]);
		dists.at(4).at(ik) = F[3][k] + F[3+6][k];
		dists.at(5).at(ik) = F[4][k] + F[4+6][k];
		dists.at(6).at(ik) = F[5][k] + F[5+6][k];
		dists.at(7).at(ik) = F[0][k];
	}

	// create the output file name to be the same as the input file
	// but with .tex instead of .dat
    // to handle windows/linux difference between specifying
	// './' or '.\', we do this:
	string::size_type dot_pos = filename.rfind('.');
	string::size_type init_pos = 0;
	if (filename[0] == '.' && (filename[1] == '/' || filename[1] == '\\'))
	{
		init_pos = 2;
		dot_pos -= 2;
	}

	string basename = filename.substr(init_pos, dot_pos);
	string name = basename + ".tex";
	
	DGLAPSolver::OutputLatexTable(dists, name, DGLAPSolver::SCIENTIFIC);
	// after making the file, copy the pdf back to the main directory
	filesystem::path pdfpath{"latex"};
	pdfpath /= basename + ".pdf";
	if (!filesystem::exists(pdfpath))
	{
		println("[ERROR] read.cpp: output pdf file \"{}\" doesnt exist...", pdfpath.string());
		exit(EXIT_FAILURE);
	}
	filesystem::copy(pdfpath, filesystem::current_path());
}

