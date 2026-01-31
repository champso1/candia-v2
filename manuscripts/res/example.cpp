#include <Candia-v2/Candia.hpp>
using namespace Candia2;

#include <vector>
#include <memory>

int main()
{
  const double num_grid_points = 1000;
  const double order = 3;
  const double Qf = 100;
  const double kr = 1.0;
  const double iterations = 10;
  const double trunc_idx = 15;

  std::vector<double> xtab{
    1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
  Grid grid(xtab, num_grid_points, 3);

  std::unique_ptr<Distribution> dist = 
    std::make_unique<LesHouchesDistribution>();
  AlphaS alphas(order, dist->Q0(), Qf, dist->alpha0(), kr);
  alphas.setVFNS(dist->masses(), dist->nfi());
  // alphas.setFFNS(4);

  DGLAPSolver solver(
    order, grid, alphas,
    Qf, iterations, trunc_idx,
    *dist, kr);
  auto dists = solver.evolve();
}