#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iostream>
#include <chrono>
using namespace std;

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"
using namespace Candia2;

int main()
{
	const double num_grid_points = 501;
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab, num_grid_points);

	vector<double> sin_grid(num_grid_points), res(num_grid_points), res_gsl(num_grid_points);
	ranges::transform(
		grid.GridPoints(), sin_grid.begin(),
		[](double x) -> double {
			return std::sin(x);
		});
	shared_ptr<Expression> P = make_shared<P3nsm>();

	chrono::high_resolution_clock clock;
	chrono::time_point<chrono::high_resolution_clock> t0, tf;
	t0 = chrono::high_resolution_clock::now();
	for (uint k=0; k<grid.Size()-1; ++k)
	{
		res[k] = grid.Convolution(sin_grid, P, k);
	}
	tf = chrono::high_resolution_clock::now();
	std::chrono::duration<double, micro> dt = tf - t0;
	std::cout << "Ordinary took " << dt.count() << " microseconds\n";

	t0 = chrono::high_resolution_clock::now();
	for (uint k=0; k<grid.Size()-1; ++k)
	{
		res[k] = grid.ConvolutionGSL(sin_grid, P, k);
	}
	tf = chrono::high_resolution_clock::now();
    dt = tf - t0;
	std::cout << "GSL took " << dt.count() << " microseconds\n";

	// ranges::for_each(res, [](double x){ cout << x << ' '; }); cout << endl;
	// ranges::for_each(res_gsl, [](double x){ cout << x << ' '; }); cout << endl;
}
