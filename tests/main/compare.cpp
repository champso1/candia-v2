#include "Candia-v2/Common.hpp"

#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Candia.hpp"
using namespace Candia2;

constexpr static char const* MSHT_DATA_FILE{"VFNS-MSHT-FHMRUVV.dat"};

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		cerr << "[ERROR] Must provide a file to read from.\n";
		exit(EXIT_FAILURE);
	}

	fs::path candia_data_path{argv[1]};
	if (!fs::exists(candia_data_path))
	{
		cerr << "[ERROR] Failed to find file " << quoted(candia_data_path.string()) << '\n';
		exit(EXIT_FAILURE);
	}

	fs::path msht_data_path{MSHT_DATA_FILE};
	if (!fs::exists(msht_data_path))
	{
		cerr << "[ERROR] Failed to find existing MSHT data in file " << quoted(MSHT_DATA_FILE) << '\n';
		exit(EXIT_FAILURE);
	}

	// fill vectors with zeros
	// 8 distributions, 11 tabulated x values each
	MultiDimVector<double, 2>::type candia_data(8, std::vector<double>(11, 0.0));
	MultiDimVector<double, 2>::type	msht_data{candia_data};
	// output vector
	MultiDimVector<double, 2>::type diffs{candia_data};
	
	// read from the candia file
	ifstream candia_stream{candia_data_path};
	candia_stream.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore the comment line

	// contains xtab and ntab lines, read those in
	vector<double> xtab{};
	vector<int> ntab{};
	double temp1{};
	int temp2;
	
	string line{};
	getline(candia_stream, line);
	istringstream iss{line};
	while (iss >> temp1)
		xtab.push_back(temp1);
	getline(candia_stream, line);
	iss = istringstream{line};
	int count = 0;
	while (iss >> temp2)
	{
	    count++;
		ntab.push_back(temp2);
	}

	// grab everything into a table F, we will parse this to match after
	MultiDimVector<double, 2>::type F(13, std::vector<double>{});
	double x{};
	while (getline(candia_stream, line))
	{
		iss = istringstream{line};
		iss >> x; // ignore this value essentially, we don't care
		for (int i=0; i<F.size(); ++i)
		{
			iss >> temp1;
			F.at(i).push_back(temp1);
		}
	}

	// candia_data now contains all the right sheisse
	// -1 on the condition since this would contain the x=1.0 but we don't want that
	for (uint ik=0; ik<ntab.size()-1; ++ik)
	{
		int k = ntab.at(ik);
		candia_data.at(0).at(ik) = F[1][k] - F[1+6][k];
		candia_data.at(1).at(ik) = F[2][k] - F[2+6][k];
		candia_data.at(2).at(ik) = F[2+6][k] - F[1+6][k];
		candia_data.at(3).at(ik) = 2.0*(F[2+6][k] + F[1+6][k]);
		candia_data.at(4).at(ik) = F[3][k] + F[3+6][k];
		candia_data.at(5).at(ik) = F[4][k] + F[4+6][k];
		candia_data.at(6).at(ik) = F[5][k] + F[5+6][k];
		candia_data.at(7).at(ik) = F[0][k];
	}

	// basically do the same but without the other shit for the msht data
	ifstream msht_stream{msht_data_path};
	if (!msht_stream)
	{
		cerr << "[ERROR] Failed to open a file stream for the MSHT data file " << quoted(MSHT_DATA_FILE) << '\n';
		exit(EXIT_FAILURE);
	}

	uint idx = 0;
	while (getline(msht_stream, line))
	{
		iss = istringstream{line};
		for (int j=0; j<msht_data.size(); ++j)
		{
			iss >> x;
			msht_data.at(j).at(idx) = x;
		}
		idx++;
	}

	auto reldiff =
		[](double candia, double msht) -> double
	{
		return abs(candia-msht)/((candia+msht)/2.0);
	};

	for (uint j=0; j<candia_data.size(); ++j)
	{
		for (uint k=0; k<candia_data.at(0).size(); ++k)
		{
			double candia = candia_data.at(j).at(k);
			double msht = msht_data.at(j).at(k);
			diffs.at(j).at(k) = reldiff(candia, msht);
		}
	}

	std::ostringstream filename_ss{};
	std::string argv1{argv[1]};
	filename_ss << "diffs-" << argv1.substr(0, argv1.find('.')) << ".tex";
	DGLAPSolver::OutputLatexTable(diffs, filename_ss.str(), DGLAPSolver::PERCENT);
}
