#include "Candia-v2/Grid.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
using namespace Candia2;

#include <fstream>
#include <array>
#include <iomanip>
#include <limits>
#include <string_view>
#include <ranges>

static void genfile();
static void compare_files(std::string_view name, std::string_view type);

static std::array ome_names{
	"AqqQNSEven",
	"AqqQNSOdd",
	"AgqQ",
	"AggQ",
	"AQqPS",
	"AQg",
	"AqqQPS",
	"AqgQ",
	"AQqPSs"
};
static std::array ome_types{
	"reg", "plus", "delta"
};

int main()
{
	genfile();
	for (auto&& name : ome_names) {
		for (auto&& type : ome_types)
			compare_files(name, type);
	}
	// compare_files("AggQ", "plus");
	// auto AggQ = ome::Candia2::AggQ;
	// std::cout << (AggQ.has_plus() ? "true" : "false") << '\n';
}

static void genfile()
{
    auto AqqQNSEven = OpMatElemN3LO(ome::Candia2::AqqQNSEven);
	auto AqqQNSOdd = OpMatElemN3LO(ome::Candia2::AqqQNSOdd);
	auto AgqQ = OpMatElemN3LO(ome::Candia2::AgqQ);
	auto AggQ = OpMatElemN3LO(ome::Candia2::AggQ);
	auto AQqPS = OpMatElemN3LO(ome::Candia2::AQqPS);
	auto AQg = OpMatElemN3LO(ome::Candia2::AQg);
	auto AqqQPS = OpMatElemN3LO(ome::Candia2::AqqQPS);
	auto AqgQ = OpMatElemN3LO(ome::Candia2::AqgQ);
	auto AQqPSs = OpMatElemN3LO(ome::Candia2::AQqPSs);
	OpMatElem::update(0, 4);

	std::array all_omes{
		std::make_pair(AqqQNSEven, "AqqQNSEven"),
		std::make_pair(AqqQNSOdd, "AqqQNSOdd"),
		std::make_pair(AgqQ, "AgqQ"),
		std::make_pair(AggQ, "AggQ"),
		std::make_pair(AQqPS, "AQqPS"),
		std::make_pair(AQg, "AQg"),
		std::make_pair(AqqQPS, "AqqQPS"),
		std::make_pair(AqgQ, "AqgQ"),
		std::make_pair(AQqPSs, "AQqPSs"),
	};
	std::array ome_types{"reg", "plus", "delta"};

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);

	namespace rv = std::views;
	auto ome_type_view =
		rv::iota(uint{0}, ome_types.size())
		| rv::transform([&](uint i){ return std::make_pair(i, ome_types[i]); });

	for (auto&& [ome,name] : all_omes) {
		for (auto&& [i,type] : ome_type_view) {
			std::ofstream outfile(std::format("data/out-new-{}-{}.dat", name, type));
			outfile << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
			for (auto [k,x] : grid.enumerate()) {
				outfile << x << ' ';
				switch (i) {
					case 0: outfile << ome.calcRegular(x) << '\n'; break;
					case 1: outfile << ome.calcPlus() << '\n'; break;
					case 2: outfile << ome.calcDelta() << '\n'; break;
				}
			}
		}
	}
}

static void compare_files(std::string_view name, std::string_view type)
{
	namespace rv = std::views;

	std::filesystem::path new_path(std::format("data/out-new-{}-{}.dat", name, type));
	std::filesystem::path old_path(std::format("data/out-old-{}-{}.dat", name, type));

	auto extract_sheisse = [](std::filesystem::path const& infile_path) {
		auto line_text = read_file(infile_path);
		auto view = std::string_view(line_text)
			| rv::split('\n')
			| rv::transform(
				[](auto&& line){
					auto line_str = std::string(line.begin(), line.end());
					std::istringstream iss_line(line_str);
					double val;
					iss_line >> val;
					iss_line >> val;
					return val;
				});
		return std::vector<double>(view.begin(), view.end());
	};

	std::vector<double> new_data = extract_sheisse(new_path);
	std::vector<double> old_data = extract_sheisse(old_path);

	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab);

	auto percent_diff = [](double x, double y) {
		double avg = std::abs(x+y)/2.0;
		return std::abs(x-y)/avg * 100.0;
	};

	std::ofstream all_output("results.txt");
	for (auto [k,x] : grid.enumerate()) {
		auto new_datapoint = new_data[k];
		auto old_datapoint = old_data[k];
		
		if (percent_diff(new_datapoint, old_datapoint) > 0.01)
			std::cout << "FAIL(" << name << "/" << type << "): x=" << x << ", new = " << new_datapoint << ", old = " << old_datapoint << '\n';
	}
}
