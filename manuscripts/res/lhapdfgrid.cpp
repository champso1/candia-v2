#include "Candia-v2/Candia.hpp"
#include "Candia-v2/LHAPDFGrid.hpp"
using namespace Candia2;

int main()
{
  std::vector<double> qvals{10.0, 50.0, 100.0};
  LesHouchesDistribution dist(qvals.back());
  Grid grid({
    1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9});
  uint order = 3;
  uint iterations = 10;
  uint trunc_idx = 10;
  double mur2_muf2 = 1.0;
  
  LHAPDFGrid lhapdfgrid(
    "testpdf", "infofile.in",
    dist, grid,
    order, iterations, trunc_idx, mur2_muf2);
  // lhapdfgrid.evolve(qvals, {});   // exact
  lhapdfgrid.evolveTrunc(qvals, {}); // truncated
  lhapdfgrid.write();
}
