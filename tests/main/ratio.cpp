#include <print>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <limits>
#include <ranges>
#include <sstream>
using namespace std;
using uint = unsigned;
namespace fs = filesystem;

static void usage()
{
	println("USAGE: ratio <n3lo-file> <nnlo-file>");
	println("    <n3lo-file>: path to n3lo data file");
	println("    <nnlo-file>: path to nnlo data file");
}

tuple<
	vector<vector<double>>, // all data
	vector<vector<double>>, // selected distributions 
	vector<int>, // ntab
	vector<double>> // all X values
readDatafile(fs::path path)
{
	print("Reading data from file \"{}\"... ", path.filename().string());
	
	ifstream file_stream{path};
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
	vector<vector<double>> F{};
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

	println("Done.");
	println("Assigning new distributions in accordance with the table... ");

	vector<vector<double>> dists(9, vector<double>{});
	double size = ntab.size();
	for (uint ix=0; ix<F.at(0).size(); ++ix)
	{
		dists[0].emplace_back(F[0].at(ix));
		dists[1].emplace_back(F[1].at(ix));
		dists[2].emplace_back(F[1+6].at(ix));
		dists[3].emplace_back(F[2].at(ix));
		dists[4].emplace_back(F[2+6].at(ix));
		dists[5].emplace_back(F[3].at(ix));
		dists[6].emplace_back(F[3+6].at(ix));
		dists[7].emplace_back(F[3].at(ix) + F[3+6].at(ix));
		dists[8].emplace_back(F[4].at(ix) + F[4+6].at(ix));
	}
	
	println("Done.");
	return {F, dists, ntab, X};
}

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		usage();
		exit(EXIT_FAILURE);
	}

	fs::path n3lo_filepath{argv[1]};
	if (!fs::exists(n3lo_filepath))
	{
		println("ratio.cpp: n3lo file given by \"{}\" doesn't exist.", n3lo_filepath.filename().string());
		exit(EXIT_FAILURE);
	}
	fs::path nnlo_filepath{argv[2]};
	if (!fs::exists(nnlo_filepath))
	{
		println("ratio.cpp: nnlo file given by \"{}\" doesn't exist.", nnlo_filepath.filename().string());
		exit(EXIT_FAILURE);
	}

	int dash_pos1{}, dash_pos2{};

	// grab the filename strings this way, rather than with directly from argv,
	// to strip away leading/training slashes (or backslashes, on Windows)
	string n3lo_filename{n3lo_filepath.filename().string()};
	dash_pos1 = n3lo_filename.find("-g");
	dash_pos2 = n3lo_filename.find("-", dash_pos1+1);
	int n3lo_gridpts = atoi(n3lo_filename.substr(dash_pos1+2, dash_pos2-dash_pos1-2).c_str());

	string nnlo_filename{nnlo_filepath.filename().string()};
	dash_pos1 = nnlo_filename.find("-g");
	dash_pos2 = nnlo_filename.find("-", dash_pos1+1);
	int nnlo_gridpts = atoi(nnlo_filename.substr(dash_pos1+2, dash_pos2-dash_pos1-2).c_str());

	if (nnlo_gridpts != n3lo_gridpts)
	{
		println("ratio.cpp: requires an equal number of grid points for both files.");
		exit(EXIT_FAILURE);
	}
   
	auto [n3lo_data_all, n3lo_data, n3lo_ntab, n3lo_xtab] = readDatafile(n3lo_filepath);
	auto [nnlo_data_all, nnlo_data, nnlo_ntab, nnlo_xtab] = readDatafile(nnlo_filepath);

	if (n3lo_data.size() != nnlo_data.size())
	{
		println("ratio.cpp: two resulting data vectors are not of the same size.");
		exit(EXIT_FAILURE);
	}
	

	vector<vector<double>> ratios{};
	for (auto const& [n3lo_v, nnlo_v] : ranges::views::zip(n3lo_data, nnlo_data))
	{
		vector<double> v{};
		for (auto const [n3lo, nnlo] : ranges::views::zip(n3lo_v, nnlo_v))
			v.emplace_back(n3lo/nnlo);
		ratios.emplace_back(v);
	}

	ofstream outfile{"ratios.dat"};
	for (int i=0; i<ratios.at(0).size()-1; ++i)
	{
		outfile << n3lo_xtab.at(i) << ' ';
		for (int j=0; j<ratios.size(); ++j)
		{
			outfile << ratios.at(j).at(i) << ' ';
		}
		outfile << '\n';
	}
}
