#include <print>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
using namespace std;
namespace fs = filesystem;
using uint = unsigned;

using vec_type = vector<vector<double>>;


static void usage()
{
	println("USAGE: ./plot-generation(.exe) <n3lo-datafile> <nnlo-datafile>");
	exit(EXIT_FAILURE);
}

tuple<vec_type, vector<double>> read_datafile(fs::path const& path);
void generate_ratios(vec_type const& n3lo_data, vec_type const& nnlo_data, vector<double> const& X);

int main(int argc, char *argv[])
{
	if (argc != 3)
		usage();

	fs::path n3lo_file{argv[1]};
	fs::path nnlo_file{argv[2]};

    auto [n3lo_data, n3lo_x] = read_datafile(n3lo_file);
	auto [nnlo_data, nnlo_x] = read_datafile(nnlo_file);

	if (n3lo_data.size() != nnlo_data.size())
	{
		println("Differing number of grid points in n3lo file ({}) and nnlo file ({}).",
			n3lo_data.at(0).size(), nnlo_data.at(0).size());
		exit(EXIT_FAILURE);
	}

	generate_ratios(n3lo_data, nnlo_data, n3lo_x);
}




tuple<vec_type, vector<double>> read_datafile(fs::path const& path)
{
	print("Reading data from file '{}'... ", path.filename().string());
	
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
	vec_type F{};
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
	return {F, X};
}

void generate_ratios(vec_type const& n3lo_data, vec_type const& nnlo_data, vector<double> const& X)
{
	ofstream outfile{"ratio.dat"};
    const uint N = n3lo_data.at(0).size();
	const uint J = n3lo_data.size();
	for (uint k=0; k<N; ++k)
	{
		outfile << X[k] << ' ';
		for (uint j=0; j<J; ++j)
			outfile << n3lo_data[j][k]/nnlo_data[j][k] << ' ';
		outfile << '\n';
	}
}
