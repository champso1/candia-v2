#include "Candia-v2/Common.hpp"
#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <filesystem>
#include <sstream>
#include <iterator>
using dist_type = std::vector<std::vector<double>>;

#include "util.hpp"

static void usage()
{
	log(LOG_INFO, "precision.cpp", "USAGE: ./precision <candia-file-1> <candia-file-2> <title>");
	log(LOG_INFO, "precision.cpp", "    <candia-file-1>: path to the first candia datafile");
	log(LOG_INFO, "precision.cpp", "    <candia-file-2>: path to the second datafile");
	log(LOG_INFO, "precision.cpp", "    <title>: title to give to the outputted table pdf");
	exit(EXIT_FAILURE);
}

static dist_type compute_diffs(dist_type const& candia, dist_type const& other);

int main(int argc, char* argv[])
{
	if (argc != 4)
		usage();

	static constexpr uint type = 1;

	fs::path filepath1(argv[1]), filepath2(argv[2]);
	file_exists(filepath1);
	file_exists(filepath2);
	std::string title(argv[3]);

    auto [xtab_candia1, candia_dists_raw1] = read_candia_file(filepath1, 13);
	auto [xtab_candia2, candia_dists_raw2] = read_candia_file(filepath2, 13);

	if (!std::ranges::equal(xtab_candia1, xtab_candia2))
	{
		log(LOG_ERROR, "precision.cpp", "Two candia datafile xtabs are different. {} vs. {}",
			vec_to_str(xtab_candia1), vec_to_str(xtab_candia2));
	}

	auto candia_dists1 = fix_dists(candia_dists_raw1, type);
	auto candia_dists2 = fix_dists(candia_dists_raw2, type);

	auto diffs = compute_diffs(candia_dists1, candia_dists2);

	std::string latex_filename = std::format("diffs-candia-{}", title);
	outputLatexTable(xtab_candia1, diffs, latex_filename, cols[type].get(), 2, false);
}

static dist_type compute_diffs(dist_type const& candia_data1, dist_type const& candia_data2)
{
	auto reldiff =
		[](double candia1, double candia2) -> double {
			double avg = (candia1+candia2)/2.0;
			return std::abs((candia1-candia2)/avg);
		};

	dist_type diffs{candia_data1.size(), std::vector<double>(candia_data1.at(0).size(), 0.0)};
	for (uint j=0; j<candia_data1.size(); ++j) {
		for (uint k=0; k<candia_data1.at(0).size(); ++k) {
			double candia1 = candia_data1.at(j).at(k);
			double candia2 = candia_data2.at(j).at(k);
			diffs.at(j).at(k) = reldiff(candia1, candia2);
		}
	}
	return diffs;
}

