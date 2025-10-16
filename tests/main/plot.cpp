#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <print>
#include <vector>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iterator>
using namespace std;
namespace fs = filesystem;
using uint = unsigned;

static void processFile(fs::path const& p);

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		println("[ERROR] please enter the LO, NLO, NNLO, and N3LO data files (in that order).");
		exit(EXIT_FAILURE);
	}

    vector<string> filenames(&argv[1], &argv[4]);
	vector<fs::path> filepaths(filenames.begin(), filenames.end());

	ranges::for_each(
		filenames,
		[](fs::path const& p) -> void {
		if (!fs::exists(p))
		{
			println("[ERROR] File \"{}\" does not exist.", p.string());
			exit(EXIT_FAILURE);
		}
	});	
}


static void processFile(fs::path const& p)
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
	while (getline(file_stream, line))
	{
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
	if (!fs::exists(plot_path) || !fs::is_directory(plot_path))
	{
		if (!fs::create_directory(plot_path))
		{
			println("[ERROR] Failed to create directory \"{}\".", plot_path.string());
			exit(EXIT_FAILURE);
		}
	}


	

	
	streampos dash_pos = p.string().find('-');
	string outfile_name = p.string().substr(0, dash_pos);
	plot_path /= fs::path{outfile_name};
	std::ofstream outfile{plot_path};

}
