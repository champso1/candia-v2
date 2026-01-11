#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/Common.hpp"
#include "LHAPDF/LHAGlue.h"
#include <charconv>
#include <cstring>
using namespace Candia2;

static void usage()
{
	log(LOG_INFO, "testing.cpp", "USAGE: ./testing <pdf> <Q0> <x>");
	log(LOG_INFO, "testing.cpp", "    <pdf>: name of LHAPDF set");
	log(LOG_INFO, "testing.cpp", "    <Q0>: energy scale");
	log(LOG_INFO, "testing.cpp", "    <x>: momentum fraction");
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	if (argc != 4) {
		log(LOG_ERROR_NOQUIT, "testing.cpp", "Invalid number of arguments: {}", argc);
		usage();
	}
	LHAPDF::setVerbosity(0);

	std::string setname(argv[1]);
	double Q0{}, x{};
	std::from_chars(argv[2], argv[2] + strlen(argv[2]), Q0);
	std::from_chars(argv[3], argv[3] + strlen(argv[3]), x);

	
	auto pdf = LHAPDFDistribution::make_lhapdf_pdf(setname, 0);
	pdf->
}
