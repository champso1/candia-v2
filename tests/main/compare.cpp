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
	log(LOG_INFO, "compare.cpp", "USAGE: ./compare <candia-file> <other-file> <format> <type>");
	log(LOG_INFO, "compare.cpp", "    <candia-file>: path to the candia datafile");
	log(LOG_INFO, "compare.cpp", "    <other-file>: path to the other datafile (should contain only tabulated values)");
	log(LOG_INFO, "compare.cpp", "    <format>: 0=benchmark format (0.1 -> 1.0^{{-1}}");
	log(LOG_INFO, "compare.cpp", "              1=normal format (0.1 -> 1e-1)");
	log(LOG_INFO, "compare.cpp", "    <type>: 0=all flavors independently");
	log(LOG_INFO, "compare.cpp", "    		  1=special combos from benchmark paper");
	log(LOG_INFO, "compare.cpp", "    		  2=special combos from benchmark paper with q(-)");
	exit(EXIT_FAILURE);
}

dist_type compute_diffs(dist_type const& candia, dist_type const& other);

int main(int argc, char *argv[])
{
	if (argc != 5)
		usage();

	fs::path candia_filepath{argv[1]}, other_filepath{argv[2]};
	file_exists(candia_filepath);
	file_exists(other_filepath);

	int type{}, format{};
	std::from_chars(argv[3], argv[3] + 1, format);
	std::from_chars(argv[4], argv[4] + 1, type);

	int ncols = cols[type].get().size();
	
	auto [xtab_candia, candia_dists_raw] = read_candia_file(candia_filepath, 13);
	auto [xtab_other, other_dists_raw] = read_other_file(other_filepath, 13);
	if (!std::ranges::equal(xtab_candia, xtab_other)) {
		log(LOG_ERROR_NOQUIT, "compare.cpp", "Two xtab arrays for the candia and other datafile are not equivalent:");

		std::ostringstream ss{};
		std::ranges::copy(xtab_candia, std::ostream_iterator<double>(ss, ", "));
		log(LOG_INFO, "compare.cpp", "Candia xtab: {}", ss.str());
		ss = {};
		std::ranges::copy(xtab_other, std::ostream_iterator<double>(ss, ", "));
		log(LOG_INFO, "compare.cpp", "Other xtab: {}", ss.str());
		exit(EXIT_FAILURE);
	}
	dist_type candia_dists = fix_dists(candia_dists_raw, type);
	dist_type other_dists = fix_dists(other_dists_raw, type);
	
	if (candia_dists_raw.at(0).size() != other_dists_raw.at(0).size()) {
		log(LOG_ERROR, "compare.cpp", "Data size mismatch. Candia size: {}, other size: {}",
			candia_dists_raw.at(0).size(), other_dists_raw.at(0).size());
	}
	if (candia_dists.size() != other_dists.size()) {
		log(LOG_ERROR, "compare.cpp", "Data size mismatch. Candia size: {}, other size: {}",
			candia_dists.size(), other_dists.size());
	}
	
	auto diffs = compute_diffs(candia_dists, other_dists);
    
	std::string identifier = candia_filepath.filename().string().substr(0, candia_filepath.filename().string().rfind('.'));
	std::string latex_filename = std::format("diffs-other-{}", identifier);
	outputLatexTable(xtab_candia, diffs, latex_filename, cols[type].get(), true, format == 0);
}


dist_type compute_diffs(dist_type const& candia_data, dist_type const& other_data)
{
	auto reldiff =
		[](double candia, double other) -> double {
			return std::abs((candia-other)/other);
		};

	dist_type diffs{candia_data.size(), std::vector<double>(candia_data.at(0).size(), 0.0)};
	for (uint j=0; j<candia_data.size(); ++j) {
		for (uint k=0; k<candia_data.at(0).size(); ++k) {
			double candia = candia_data.at(j).at(k);
			double other = other_data.at(j).at(k);
			diffs.at(j).at(k) = reldiff(candia, other);
		}
	}
	return diffs;
}
