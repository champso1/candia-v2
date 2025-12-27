#include <print>
#include <tuple>
#include <vector>
#include <ranges>

#include "Candia-v2/OperatorMatrixElements.hpp"
using ome_type = ome::rpd_distribution<ome::ome_as_view<double>, ome::ome_as_plus_view<double>, ome::ome_as_const_view<double>>;
using namespace Candia2;

static constexpr uint N = 1000;
static uint NF = 3;
static double LM = 0.0;

static A2ns a2ns_candia{};
static ome_type a2ns_libome = ome::AqqQNSEven;

static A2gq a2gq_candia{};
static ome_type a2gq_libome = ome::AgqQ;

static A2gg a2gg_candia{};
static ome_type a2gg_libome = ome::AggQ;

static A2hq a2hq_candia{};
static ome_type a2hq_libome = ome::AQqPS;

static A2hg a2hg_candia{};
static ome_type a2hg_libome = ome::AQg;

static void update(uint nf, double lm)
{
	OpMatElem::update(lm, nf);
	NF = nf;
}

static std::tuple<double,double,double> get_candia(Expression& ome, double x)
{
	return {ome._reg_func(x), ome._plus_func(x), ome._delta_func(x)};
}

static std::tuple<double,double,double> get_libome(ome_type& ome, double x)
{
	return {ome.has_regular() ? ome.get_regular().value()[2](LM, NF, x) : 0.0,
			ome.has_plus() ? ome.get_plus().value()[2](LM, NF, x) : 0.0,
			ome.has_delta() ? ome.get_delta().value()[2](LM, NF) : 0.0};
}

int main()
{
	update(3, 0.0);
	std::vector<double> grid_points = std::ranges::views::iota(1)
		| std::ranges::views::take(N)
		| std::ranges::views::transform([=](int i) -> double { return static_cast<double>(i)/N; })
		| std::ranges::to<std::vector<double>>();

	std::vector<std::tuple<double,double,double>> candia_results{}, libome_results{}, diffs{};
	for (double x : grid_points)
	{
		candia_results.emplace_back(get_candia(a2hg_candia, x));
		libome_results.emplace_back(get_libome(a2hg_libome, x));

	    diffs.emplace_back(
			std::get<0>(candia_results.back())/std::get<0>(libome_results.back()),
			std::get<1>(candia_results.back()) == 0.0 ? 0.0 : std::get<1>(candia_results.back())/std::get<1>(libome_results.back()),
			std::get<2>(candia_results.back()) == 0.0 ? 0.0 : std::get<2>(candia_results.back())/std::get<2>(libome_results.back()));
	}
	
    for (auto t : diffs)
		std::println("{}", t);
}
