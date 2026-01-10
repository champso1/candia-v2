#include <cstdlib>
#include <ranges>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;
using uint = unsigned;

using vec_type = std::vector<std::vector<double>>;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

static void usage()
{
	log(LOG_INFO, "plot-generation", "USAGE: ./plot-generation <n3lo-datafile> <nnlo-datafile>");
	exit(EXIT_FAILURE);
}

std::tuple<vec_type, std::vector<double>> read_datafile(fs::path const& path);
void generate_ratios(vec_type const& n3lo_data, vec_type const& nnlo_data, std::vector<double> const& X);

int main(int argc, char *argv[])
{
	if (argc != 3) {
		log(LOG_ERROR_NOQUIT, "plot-generation", "Invalid number of arguments: {}", argc-1);
		usage();
	}

	fs::path n3lo_file{argv[1]};
	fs::path nnlo_file{argv[2]};

    auto [n3lo_data, n3lo_x] = read_datafile(n3lo_file);
	auto [nnlo_data, nnlo_x] = read_datafile(nnlo_file);

	if (n3lo_data.size() != nnlo_data.size()) {
		log(LOG_ERROR, "plot-generation", "Differing number of grid points in n3lo file ({}) and nnlo file ({}).",
			n3lo_data.at(0).size(), nnlo_data.at(0).size());
		exit(EXIT_FAILURE);
	}

	if (auto res = std::ranges::mismatch(n3lo_x, nnlo_x); res.in1 != n3lo_x.end()) {
		double n3lo_val = *res.in1;
		double nnlo_val = *res.in2;
		log(LOG_ERROR, "plot-generation", "x values differ: n3lo val={}, nnlo val={}", n3lo_val, nnlo_val);
		exit(EXIT_FAILURE);
	}

	generate_ratios(n3lo_data, nnlo_data, n3lo_x);
}




std::tuple<vec_type, std::vector<double>> read_datafile(fs::path const& path)
{
	log(LOG_INFO, "plot-generation", "Reading data from file '{}'... ", path.filename().string());
	
	std::ifstream file_stream{path};
	file_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // ignore the comment line

	std::vector<double> xtab{};
	std::vector<int> ntab{};
	
	double temp1{};
	int temp2{};
	std::string line{};

	// read in xtab array
	std::getline(file_stream, line);
	std::istringstream iss{line};
	while (iss >> temp1)
		xtab.push_back(temp1);

	// read in ntab array
	getline(file_stream, line);
	iss = std::istringstream{line};
	while (iss >> temp2)
		ntab.push_back(temp2);

	// read in rest of data points
	std::vector<double> X{};
	vec_type F{};
	F.resize(13);
	while (getline(file_stream, line)) {
		iss = std::istringstream{line};
		iss >> temp1;
		if (temp1 > 0.95)
			break;
		X.push_back(temp1);
		for (int i=0; i<F.size(); ++i) {
			iss >> temp1;
			F.at(i).push_back(temp1);
		}
	}

	vec_type dists(6, std::vector<double>(F.at(0).size(), 0.0));
	for (uint k=0; k<F.at(0).size(); ++k) {
		dists[0][k] = F[5][k] + F[6+5][k];
		dists[1][k] = F[4][k] + F[6+4][k];
		dists[2][k] = F[3][k] + F[6+3][k];
		dists[3][k] = F[1][k] - F[6+1][k];
		dists[4][k] = F[2][k] - F[6+2][k];
		dists[5][k] = F[0][k];
	}

	return {dists, X};
}

void generate_ratios(vec_type const& n3lo_data, vec_type const& nnlo_data, std::vector<double> const& X)
{
	std::ofstream outfile{"ratio.dat"};
    const uint N = n3lo_data.at(0).size();
	const uint J = n3lo_data.size();
	for (uint k=0; k<N; ++k) {
		outfile << X[k] << ' ';
		for (uint j=0; j<J; ++j)
			outfile << n3lo_data[j][k]/nnlo_data[j][k] << ' ';
		outfile << '\n';
	}
}
