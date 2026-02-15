#include "Candia-v2/Common.hpp"
#include "Candia-v2/SplittingFn.hpp"
using namespace Candia2;

#include <ranges>
#include <memory>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <functional>

struct splitfunc_type_t final
{
	using value_type = std::unique_ptr<SplittingFunction>;
	value_type splitfunc{nullptr};
	std::string_view title{"P(x)"};
};

template <typename TSplitFunc>
static splitfunc_type_t make_splitfunc_type(std::string_view title)
{
	return splitfunc_type_t{
		.splitfunc = std::make_unique<TSplitFunc>(),
		.title = title
	};
}

int main()
{
	std::vector<splitfunc_type_t> splitfuncs{};
    splitfuncs.emplace_back(make_splitfunc_type<P0ns>("p0ns"));
	splitfuncs.emplace_back(make_splitfunc_type<P0qq>("p0qq"));
	splitfuncs.emplace_back(make_splitfunc_type<P0qg>("p0qg"));
	splitfuncs.emplace_back(make_splitfunc_type<P0gq>("p0gq"));
	splitfuncs.emplace_back(make_splitfunc_type<P0gg>("p0gg"));
	splitfuncs.emplace_back(make_splitfunc_type<P1nsm>("p1nsm"));
	splitfuncs.emplace_back(make_splitfunc_type<P1nsp>("p1nsp"));
	splitfuncs.emplace_back(make_splitfunc_type<P1qq>("p1qq"));
	splitfuncs.emplace_back(make_splitfunc_type<P1qg>("p1qg"));
	splitfuncs.emplace_back(make_splitfunc_type<P1gq>("p1gq"));
	splitfuncs.emplace_back(make_splitfunc_type<P1gg>("p1gg"));
	splitfuncs.emplace_back(make_splitfunc_type<P2nsm>("p2nsm"));
	splitfuncs.emplace_back(make_splitfunc_type<P2nsp>("p2nsp"));
	splitfuncs.emplace_back(make_splitfunc_type<P2nsv>("p2nsv"));
	splitfuncs.emplace_back(make_splitfunc_type<P2qq>("p2qq"));
	splitfuncs.emplace_back(make_splitfunc_type<P2qg>("p2qg"));
	splitfuncs.emplace_back(make_splitfunc_type<P2gq>("p2gq"));
	splitfuncs.emplace_back(make_splitfunc_type<P2gg>("p2gg"));
	splitfuncs.emplace_back(make_splitfunc_type<P3nsm>("p3nsm"));
	splitfuncs.emplace_back(make_splitfunc_type<P3nsp>("p3nsp"));
	splitfuncs.emplace_back(make_splitfunc_type<P3nsv>("p3nsv"));
	splitfuncs.emplace_back(make_splitfunc_type<P3qq>("p3qq"));
	splitfuncs.emplace_back(make_splitfunc_type<P3qg>("p3qg"));
	splitfuncs.emplace_back(make_splitfunc_type<P3gq>("p3gq"));
	splitfuncs.emplace_back(make_splitfunc_type<P3gg>("p3gg"));

	std::filesystem::path template_gplt_file_path("plot-scripts/all-split-funcs.gplt");
	if (!std::filesystem::exists(template_gplt_file_path))
		log(LOG_ERROR, "generate-splitfunc-plots.cpp", "Missing template file: {}", template_gplt_file_path.string());
	std::ifstream template_gplt_file(template_gplt_file_path);
	std::ostringstream template_gplt_ss;
	std::copy(
		std::istreambuf_iterator<char>(template_gplt_file),
		std::istreambuf_iterator<char>{},
		std::ostreambuf_iterator<char>(template_gplt_ss));
	std::string template_gplt{std::move(template_gplt_ss.str())};

    auto pipe_deleter = [](FILE* f){ pclose(f); };
    using pipe_deleter_type = decltype(pipe_deleter);
	using pipe_type = std::unique_ptr<FILE, pipe_deleter_type>;
	pipe_type gnuplot = pipe_type(popen("gnuplot", "w"), pipe_deleter);
	fprintf(gnuplot.get(), "%s", template_gplt.c_str());

    uint num_points = 10000;
	double logmin = std::log10(1e-7);
	double logmax = std::log10(1.0);
	double linmin = 0.75;
	double linmax = 1.0;
	auto points_log =
		std::views::iota(0)
		| std::views::take(num_points)
		| std::views::transform([&](int i){ return std::pow(10.0, (logmin + (logmax-logmin)*static_cast<double>(i)/static_cast<double>(num_points))); });
	auto points_lin =
		std::views::iota(0)
		| std::views::take(num_points)
		| std::views::transform([&](int i){ return linmin + (linmax-linmin)*static_cast<double>(i)/static_cast<double>(num_points); });
	std::filesystem::path datafile_dir_path("split-func-data");
	if (!std::filesystem::exists(datafile_dir_path)) {
		if (!std::filesystem::create_directory(datafile_dir_path))
			log(LOG_ERROR, "generate-splitfunc-plots.cpp", "Failed to create datafile directory: {}", datafile_dir_path.string());
	}

	auto calc = [&](auto&& points, splitfunc_type_t& splitfunc) {
		log(LOG_INFO, "generate-splitfunc-plots.cpp", "calc: {}", splitfunc.title);
		std::string datafile_name = std::format("{}.dat", splitfunc.title);
		std::filesystem::path datafile_path = datafile_dir_path/datafile_name;
		std::ofstream datafile(datafile_path);

		for (double x : points) {
			datafile << x << ' '
					 << splitfunc.splitfunc->calcRegular(x) << ' '
					 << splitfunc.splitfunc->calcPlus(x) << ' '
					 << splitfunc.splitfunc->calcDelta(x) << '\n';
		}
	};

	std::filesystem::path output_dir_path("split-func-out");
	if (!std::filesystem::exists(output_dir_path)) {
		if (!std::filesystem::create_directory(output_dir_path))
			log(LOG_ERROR, "generate-splitfunc-plots.cpp", "Failed to create output directory: {}", output_dir_path.string());
	}
	std::string_view set_output_cmd_template = "set output 'split-func-out/%NAME%-%TYPE%.png'\n";
	std::string_view plot_cmd_template = "plot 'split-func-data/%NAME%.dat' using 1:%IDX% with lines lw 2\n";
	std::string_view title_cmd_template = "set title '%NAME% %TYPE%'\n";
	std::vector<std::pair<std::string_view, std::string_view>> splitfunc_types{
		{"reg", "2"},
		{"plus", "3"},
		{"delta", "4"}
	};
	auto plot = [&](splitfunc_type_t& splitfunc) {
		log(LOG_INFO, "generate-splitfunc-plots.cpp", "plot: {}", splitfunc.title);
	    auto replace = [&](
			std::string& s, std::string_view repl, std::string_view val)
		{
			auto it = s.find(repl);
			while (it != std::string::npos) {
				s.replace(it, repl.size(), val);
				it = s.find(repl);
			}
		};

		for (auto [type, idx] : splitfunc_types) {
			std::string set_output_cmd(set_output_cmd_template);
			std::string plot_cmd(plot_cmd_template);
			std::string title_cmd(title_cmd_template);
			
			replace(set_output_cmd, "%NAME%", splitfunc.title);
			replace(set_output_cmd, "%TYPE%", type);
		    replace(plot_cmd, "%NAME%", splitfunc.title);
			replace(plot_cmd, "%IDX%", idx);
			replace(title_cmd, "%NAME%", splitfunc.title);
			replace(title_cmd, "%TYPE%", type);

			fprintf(gnuplot.get(), "%s", title_cmd.c_str());
			fprintf(gnuplot.get(), "%s", set_output_cmd.c_str());
			fprintf(gnuplot.get(), "%s", plot_cmd.c_str());

			log(LOG_INFO, "generate-splitfunc-plots.cpp", "  - Plot command ({}): {}", type, set_output_cmd);
			
			fprintf(gnuplot.get(), "clear\n");
		}
	};

	using namespace std::placeholders;
	auto calc_log = std::bind(calc, points_log, _1);
	auto calc_lin = std::bind(calc, points_lin, _1);

	std::ranges::for_each(splitfuncs, calc_lin);
	std::ranges::for_each(splitfuncs, plot);

	fprintf(gnuplot.get(), "exit\n");
}
