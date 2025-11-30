#include <memory>
#include <print>
#include <vector>
using uint = unsigned;

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/Candia.hpp"

int main()
{
	const uint num_grid_points = 100;
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Candia2::Grid grid(xtab, num_grid_points);
	
	const double Qf = 100.0;
	Candia2::DGLAPSolver solver(2, grid, Qf, 5, 5, std::make_unique<Candia2::LesHouchesDistribution>(), 1.0);
	auto dists = solver.Evolve();
	
	for (auto const [i, x] : std::ranges::views::enumerate(grid.Points()))
		std::println("{:.3e}   {:.3e}  {:.3e}  {:.3e}", x, dists[0][i], dists[1][i], dists[2][i]);
}

