#include <cstdlib>
#include <charconv>
#include <string>
#include <vector>
#include <filesystem>
using uint = unsigned;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

#include "util.hpp"
using value_type = double;

static void usage()
{
	log(LOG_INFO, "read.cpp", "USAGE: ./read <data-file> <origin> <format> <type>");
	log(LOG_INFO, "read.cpp", "    <origin>: 0=candia");
	log(LOG_INFO, "read.cpp", "              1=other (only tabulated x-values)");
	log(LOG_INFO, "read.cpp", "    <format>: 0=benchmark format (0.1 -> 1.0^{{-1}}");
	log(LOG_INFO, "read.cpp", "              1=normal format (0.1 -> 1e-1)");
	print_compare_types("read.cpp");
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	if (argc != 5 && argc != 2)
		usage();

	fs::path datafile_path(argv[1]);
	file_exists(datafile_path);

	int origin{}, type{}, format{};
	if (argc != 2) {
		std::from_chars(argv[2], argv[2] + 1, origin);
		std::from_chars(argv[3], argv[3] + 1, format);
		std::from_chars(argv[4], argv[4] + 1, type);
	} else {
		origin = 0;
		format = 1;
		type = 2;
	}

	auto [xtab, dists_raw] =
		(origin == 0) ?
		read_candia_file<value_type>(datafile_path, 37)
		: read_other_file<value_type>(datafile_path, 37);
    auto dists = fix_dists(dists_raw, type);

	std::string basename = datafile_path.filename().string().substr(0, datafile_path.filename().string().rfind('.'));	
	outputLatexTable(xtab, dists, basename, cols[type].get(), 1, format == 0);
}
