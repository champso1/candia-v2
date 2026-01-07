#include <cstdlib>
#include <string>
#include <vector>
#include <filesystem>
using uint = unsigned;
using dist_type = std::vector<std::vector<double>>;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

#include "util.hpp"


static std::vector<double> XTAB{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};

static void usage()
{
	log(LOG_INFO, "read.cpp", "USAGE: ./read <candia-file> <origin> <type>");
	log(LOG_INFO, "read.cpp", "    <origin>: 0=candia");
	log(LOG_INFO, "read.cpp", "              1=other (only tabulated x-values)");
	log(LOG_INFO, "read.cpp", "    <type>: 0=all flavors independently");
	log(LOG_INFO, "read.cpp", "            1=special combos from benchmark paper");
	log(LOG_INFO, "read.cpp", "            2=special combos from benchmark paper, with q(-)");
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	if (argc != 4)
		usage();

	fs::path datafile_path(argv[1]);
	file_exists(datafile_path);

	int origin{}, type{};
	std::from_chars(argv[2], argv[2] + 1, origin);
	std::from_chars(argv[3], argv[3] + 1, type);
	
	dist_type dists_raw = origin == 0 ? read_candia_file(datafile_path, 13) : read_other_file(datafile_path, 13);
	dist_type dists = fix_dists(dists_raw, type);

	std::string basename = datafile_path.filename().string().substr(0, datafile_path.filename().string().rfind('.'));	
	outputLatexTable(XTAB, dists, basename, cols[type].get(), false);
}
