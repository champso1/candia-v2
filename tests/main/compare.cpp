#include "util.hpp"

static void usage()
{
	log(LOG_INFO, "compare.cpp", "USAGE: ./compare <candia-file> <other-file> <format> <type> <diff-type>");
	log(LOG_INFO, "compare.cpp", "    <candia-file>: path to the candia datafile");
	log(LOG_INFO, "compare.cpp", "");
	log(LOG_INFO, "compare.cpp", "    <other-file>: path to the other datafile (should contain only tabulated values)");
	log(LOG_INFO, "compare.cpp", "");
	log(LOG_INFO, "compare.cpp", "    <format>:    0: benchmark format (0.1 -> 1.0^{{-1}}");
	log(LOG_INFO, "compare.cpp", "                 1: normal format (0.1 -> 1e-1)");
	log(LOG_INFO, "compare.cpp", "");
	print_compare_types("compare.cpp");
	log(LOG_INFO, "compare.cpp", "");
	log(LOG_INFO, "compare.cpp", "    <diff-type>: 0: percent error (i.e treating the other file as the base)");
	log(LOG_INFO, "compare.cpp", "                 1: percent difference (i.e. treating neither file as the base)");
	exit(EXIT_FAILURE);
}

static dist_type compute_diffs(
	dist_type const& candia,
	dist_type const& other,
	int diff_type);

int main(int argc, char *argv[])
{
	if (argc != 6 && argc != 3)
		usage();

	fs::path candia_filepath{argv[1]}, other_filepath{argv[2]};
	file_exists(candia_filepath);
	file_exists(other_filepath);

	int format{}, type{}, diff_type{};
	if (argc == 3) {
		format = 1;
		type = 2;
		diff_type = 1;
	} else {
		std::from_chars(argv[3], argv[3] + 1, format);
		std::from_chars(argv[4], argv[4] + 1, type);
		std::from_chars(argv[5], argv[5] + 1, diff_type);
	}

	auto read_candia_file_result = read_candia_file(candia_filepath, 37);
	auto const& xtab_candia = read_candia_file_result.xtab;
	auto const& candia_dists_raw = read_candia_file_result.dists_ntabbed;
	[[maybe_unused]] auto const& grid_points = read_candia_file_result.grid_points;
	auto [xtab_other, other_dists_raw] = read_other_file(other_filepath, 37);
	if (!std::ranges::equal(xtab_candia, xtab_other)) {
		log(LOG_ERROR_NOQUIT, "compare.cpp", "Two xtab arrays for the candia and other datafile are not equivalent:");

		log(LOG_INFO, "compare.cpp", "Candia xtab: {}", vec_to_str(xtab_candia));
		log(LOG_INFO, "compare.cpp", "Other xtab: {}", vec_to_str(xtab_other));
		exit(EXIT_FAILURE);
	}
	auto candia_dists = fix_dists(candia_dists_raw, type);
	auto other_dists = fix_dists(other_dists_raw, type);
	
	if (candia_dists_raw.at(0).size() != other_dists_raw.at(0).size()) {
		log(LOG_ERROR, "compare.cpp", "Data size mismatch. Candia size: {}, other size: {}",
			candia_dists_raw.at(0).size(), other_dists_raw.at(0).size());
	}
	if (candia_dists.size() != other_dists.size()) {
		log(LOG_ERROR, "compare.cpp", "Data size mismatch. Candia size: {}, other size: {}",
			candia_dists.size(), other_dists.size());
	}
	
	auto diffs = compute_diffs(candia_dists, other_dists, diff_type);
    
	std::string identifier = candia_filepath.filename().string().substr(0, candia_filepath.filename().string().rfind('.'));
	std::string latex_filename = std::format("diffs-other-{}", identifier);
	outputLatexTable(xtab_candia, diffs, latex_filename, cols[type].get(), 0, format == 0);
}


static dist_type compute_diffs(
	dist_type const& candia_data,
	dist_type const& other_data,
	int diff_type)
{
	auto reldiff =
		[&](double candia, double other) -> double {
			double base = diff_type == 0 ? other : (candia+other)/2.0;
			return std::abs((candia-other)/base)*100.0;
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
