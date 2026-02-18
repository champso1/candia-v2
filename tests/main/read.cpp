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
	log(LOG_INFO, "read.cpp", "USAGE: ./read <candia-file> <origin> <format> <type>");
	log(LOG_INFO, "read.cpp", "    <origin>: 0=candia");
	log(LOG_INFO, "read.cpp", "              1=other (only tabulated x-values)");
	log(LOG_INFO, "read.cpp", "    <format>: 0=benchmark format (0.1 -> 1.0^{{-1}}");
	log(LOG_INFO, "read.cpp", "              1=normal format (0.1 -> 1e-1)");
	log(LOG_INFO, "read.cpp", "    <type>: 0=all flavors independently");
	log(LOG_INFO, "read.cpp", "            1=special combos from benchmark paper");
	log(LOG_INFO, "read.cpp", "            2=special combos from benchmark paper, with q(-)");
	log(LOG_INFO, "read.cpp", "            3=mixture of special singlet and non-singlet");
	log(LOG_INFO, "read.cpp", "            4=ffns");
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	if (argc != 5)
		usage();

	fs::path datafile_path(argv[1]);
	file_exists(datafile_path);

	int origin{}, type{}, format{};
	std::from_chars(argv[2], argv[2] + 1, origin);
	std::from_chars(argv[3], argv[3] + 1, format);
	std::from_chars(argv[4], argv[4] + 1, type);

	auto [xtab, dists_raw] =
		(origin == 0) ?
		read_candia_file<value_type>(datafile_path, 13)
		: read_other_file<value_type>(datafile_path, 13);
	dist_type dists = fix_dists(dists_raw, type);

	std::string basename = datafile_path.filename().string().substr(0, datafile_path.filename().string().rfind('.'));	
	outputLatexTable(xtab, dists, basename, cols[type].get(), 1, format == 0);
}
