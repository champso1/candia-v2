#include "Candia-v2/Common.hpp"
#include <filesystem>
#include <iomanip>
using namespace Candia2;

#include "algorithm"

#include "util.hpp"

static void usage()
{
	std::cout << "[ERROR] ratio.cpp: Invalid arguments.\n";
	std::cout << "Usage:\n";
	std::cout << "-------------------------------------------------------\n";
	std::cout << "./ratio(.exe) <file-base> <file1> <title1> [<fileN> <titleN> ...]\n";
	std::cout << "    <file-base>: the base file the others are ratios to.\n";
	std::cout << "    <file1>:     the first non-base file to perform ratios to\n";
	std::cout << "    <title1>:    the title of the output datafile (no extension) for data1/data-base\n";
	std::cout << "    <[fileN]>:   additional files\n";
	std::cout << "    <[titleN]>:  additional titles\n";
	std::cout << "-------------------------------------------------------\n\n";
	std::cout << "NOTE: will emit a datafile for each non-base file provided" << std::endl;
	throw std::runtime_error("invalid cli arguments");
}


int main(int argc, char *argv[])
{
	getLogOptions().verbosity = LOG_INFO;
	if (argc < 4 || (argc%2 != 0))
		usage();

	fs::path filepath_base(argv[1]);
	if (!fs::exists(filepath_base))
		log(LOG_ERROR, "ratio.cpp", "base file doesn't exist ({})", filepath_base.string());
	ReadCandiaFileResult basefile_data = read_candia_file(filepath_base, 37);
	
	
	std::vector<fs::path> filepaths{};
	std::vector<std::string> titles{};
	std::vector<ReadCandiaFileResult> addfiles_results{};
	for (uint i=4; i<argc; i+=2) {
		filepaths.emplace_back(argv[i]);
		if (!fs::exists(filepaths.back()))
			log(LOG_ERROR, "ratio.cpp", "additional file doesn't exist ({})", filepaths.back().string());
		titles.emplace_back(argv[i+1]);
		addfiles_results.emplace_back(read_candia_file(filepaths.back(), 37));
	}

	if (
		std::ranges::any_of(addfiles_results, [&](ReadCandiaFileResult const& r){
			if (!std::ranges::equal(r.xtab, basefile_data.xtab))
				return true;
			if (r.dists_raw.front().size() != basefile_data.dists_raw.front().size())
				return true;
			// maybe more idk
			return false;
		})
	) {
		log(LOG_ERROR, "ratio.cpp", "mismatch of some kind in the datafiles");
	}

	fs::path outdir = fs::current_path()/fs::path("ratios");
	if (!fs::exists(outdir) || !fs::is_directory(outdir)) {
		if (!fs::create_directory(outdir))
			log(LOG_ERROR, "ratio.cpp", "failed to create output dir (or doesn't exist): {}", outdir.string());
	}

	{
		fs::path datafile = outdir/fs::path("base.dat");
		std::ofstream datastream(datafile);
		
		datastream << std::scientific;

		dist_type ratios = basefile_data.dists_raw; // copy
		for (uint j=0; j<ratios.size(); ++j) {
			for (uint k=0; k<ratios.front().size(); ++k)
				ratios[j][k] /= basefile_data.dists_raw[j][k];
		}

		for (uint k=0; k<ratios.front().size(); ++k) {
			datastream << std::setprecision(1) << basefile_data.xtab[k] << ' ';
			datastream << std::setprecision(std::numeric_limits<double>::max_digits10);
			for (uint j=0; j<ratios.size(); ++j)
				datastream << ratios[j][k] << ' ';
			datastream << std::endl;
		}
	}

	auto enum_view =
		std::views::iota(uint{0}, titles.size())
		| std::views::transform([&](uint i){ return std::make_pair(addfiles_results[i].dists_raw, titles[i]); });
	for (auto&& _t : enum_view) {
		dist_type const& dist = _t.first;
		std::string const& title = _t.second;
		
		std::string outfile = std::format("{}.dat", title);
		fs::path datafile = outdir/fs::path(outfile);
		std::ofstream datastream(datafile);
		datastream << std::scientific;

		dist_type ratios = dist; // copy
		for (uint j=0; j<ratios.size(); ++j) {
			for (uint k=0; k<ratios.front().size(); ++k)
				ratios[j][k] /= basefile_data.dists_raw[j][k];
		}

		for (uint k=0; k<ratios.front().size(); ++k) {
			datastream << std::setprecision(1) << basefile_data.xtab[k] << ' ';
			datastream << std::setprecision(std::numeric_limits<double>::max_digits10);
			for (uint j=0; j<ratios.size(); ++j)
				datastream << ratios[j][k] << ' ';
			datastream << std::endl;
		}
	}
}
