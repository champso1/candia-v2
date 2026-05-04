#include <algorithm>
#include <filesystem>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

[[maybe_unused]]
static void process_file(fs::path const& p);

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		log(LOG_ERROR, "plot.cpp", "Please enter the LO, NLO, NNLO, and N3LO data files (in that order).");
		exit(EXIT_FAILURE);
	}

    vector<string> filenames(&argv[1], &argv[4]);
	vector<fs::path> filepaths(filenames.begin(), filenames.end());

	ranges::for_each(
		filenames,
		[](fs::path const& p) -> void {
		if (!fs::exists(p))
			log(LOG_ERROR, "plot.cpp", "File '{}' does not exist.", p.string());
	});	
}


static void process_file(fs::path const& p)
{
	vector<double> X{};
	vector<vector<double>> F{};
	ifstream file_stream{p};
	// ignore first three lines: comment, xtab, and ntab
	for (uint i=0; i<3; ++i)
		file_stream.ignore(numeric_limits<streamsize>::max(), '\n');

	string line{};
	istringstream iss{};
	double x{};
	while (getline(file_stream, line)) {
		iss = istringstream{line};
		iss >> x; // X value
		X.emplace_back(x);

		// rest of the elements
		F.emplace_back(
			istream_iterator<double>{iss},
			istream_iterator<double>{}
		);
	}

	// make a temporary plots directly to save data file in
	fs::path plot_path{"plot"};
	if (!fs::exists(plot_path) || !fs::is_directory(plot_path)) {
		if (!fs::create_directory(plot_path))
			log(LOG_ERROR, "plot.cpp", "Failed to create directory '{}'.", plot_path.string());
	}
	
	streampos dash_pos = p.string().find('-');
	string outfile_name = p.string().substr(0, dash_pos);
	plot_path /= fs::path{outfile_name};
	std::ofstream outfile{plot_path};
}
