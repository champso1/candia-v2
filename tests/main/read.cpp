#include "Candia-v2/Common.hpp"
using namespace Candia2;

#include "util.hpp"

static void usage()
{
	log(LOG_INFO, "read.cpp", "USAGE: ./read <data-file> <origin> <format> <type>");
	log(LOG_INFO, "read.cpp", "    <origin>: 0=candia");
	log(LOG_INFO, "read.cpp", "              1=other (only tabulated x-values)");
	log(LOG_INFO, "read.cpp", "    <format>: 0=benchmark format (0.1 -> 1.0^{{-1}}");
	log(LOG_INFO, "read.cpp", "              1=normal format (0.1 -> 1e-1)");
	print_compare_types("read.cpp");
	throw std::runtime_error("invalid cli arguments");
}

int main(int argc, char *argv[])
{
	getLogOptions().verbosity = LOG_INFO;
	if (argc != 5 && argc != 2)
		usage();

	fs::path datafile_path(argv[1]);
	file_exists(datafile_path);

	int origin{}, format{}, type{};
	if (argc == 2) {
		origin = 0;
		format = 1;
		type = 1;
	} else {
		std::from_chars(argv[2], argv[2] + 1, origin);
		std::from_chars(argv[3], argv[3] + 1, format);
		std::from_chars(argv[4], argv[4] + 1, type);
	}

	xtab_type xtab{};
	dist_type dists_raw{};
	if (origin == 0 && type != 5) {
		auto read_candia_file_result = read_candia_file(datafile_path, 37);
		xtab = read_candia_file_result.xtab;
		dists_raw = read_candia_file_result.dists_ntabbed;
	} else if (origin != 0 && type != 5){
		auto [xtab, dists] = read_other_file(datafile_path, 37);
	} else if (type == 5) {
		auto num_dists = static_cast<uint>(QEDPartonIndices::COUNT);
		auto read_candia_file_result = read_candia_file(datafile_path, num_dists);
		xtab = read_candia_file_result.xtab;
		dists_raw = read_candia_file_result.dists_ntabbed;
	}

	auto dists = fix_dists(dists_raw, type);
	
	std::string basename = datafile_path.filename().string().substr(0, datafile_path.filename().string().rfind('.'));	
	outputLatexTable(xtab, dists, basename, cols[type].get(), 1, format == 0);
}
