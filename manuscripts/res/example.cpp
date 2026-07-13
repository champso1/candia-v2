#include <Candia-v2/Candia.hpp>
using namespace Candia2;

#include <vector>

int main()
{
  const double order = 3;
  const double iterations = 10;
  const double trunc_idx = 10;
  const double Qf = 100;
  const double mur2_muf2 = 1.0;

  std::vector<double> xtab{
    1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
  Grid grid(xtab);

  LesHouchesDistribution dist(Qf);
  AlphaS alphas(
    order,
    dist.Q0(), dist.Qf(),
    dist.alpha0(),
    mur2_muf2);
  alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
  // alphas.setFFNS(4);

  DGLAPSolver solver(
    order, grid, alphas,
    Qf, iterations, trunc_idx,
    dist, mur2_muf2);
  std::vector<ArrayGrid> dists = solver.evolve();
  // do stuff with the distributions
}
