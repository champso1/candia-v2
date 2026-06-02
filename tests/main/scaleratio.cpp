#include "Candia-v2/Common.hpp"
using namespace Candia2;

#include "util.hpp"

static void usage()
{
	log(LOG_INFO, "scaleratio.cpp", "USAGE: ./scaleratio <data-file-mu1> <data-file-mu0.5> <data-file-mu2> <format> <type> <title>");
	log(LOG_INFO, "scaleratio.cpp", "    <data-file-mu[1|0.5|2]>: datafile corresponding to mur_muf=[1|0.5|2]");
	log(LOG_INFO, "scaleratio.cpp", "    <format>: 0=benchmark format (0.1 -> 1.0^{{-1}}");
	log(LOG_INFO, "scaleratio.cpp", "              1=normal format (0.1 -> 1e-1)");
	print_compare_types("scaleratio.cpp");
	log(LOG_INFO, "scaleratio.cpp", "    <title>: title for single resulting table");
	throw std::runtime_error("invalid cli arguments");
}

int main(int argc, char *argv[])
{
	getLogOptions().verbosity = LOG_INFO;
	if (argc != 7 && argc != 5)
		usage();

	std::string datafile_mu1{argv[1]};
	std::string datafile_mu05{argv[2]};
	std::string datafile_mu2{argv[3]};
	int format{};
	int type{};
	std::string title{argv[6]};

	if (argc == 4) {
		format = 1;
		type = 2;
	} else {
		std::from_chars(argv[4], argv[4] + 1, format);
		std::from_chars(argv[5], argv[5] + 1, type);
	}

	xtab_type xtab{};
	dist_type dists_raw{};
	auto mu1_result = read_candia_file(datafile_mu1, 37);
	auto mu05_result = read_candia_file(datafile_mu05, 37);
	auto mu2_result = read_candia_file(datafile_mu2, 37);
	
	auto mu1_dists = fix_dists(mu1_result.dists_ntabbed, type);
	auto mu05_dists = fix_dists(mu05_result.dists_ntabbed, type);
	auto mu2_dists = fix_dists(mu2_result.dists_ntabbed, type);

	std::array<std::reference_wrapper<dist_type>, 3> arrs{
		std::ref(mu1_dists),
		std::ref(mu05_dists),
		std::ref(mu2_dists)
	};
	
	outputLatexTableScaleRatio(mu1_result.xtab, arrs, title, cols[type].get(), 1, format == 0);
}

