// numerical-precision.cpp
//
// Measures the numerical error of a single convolution (E (x) A)(x) as Candia computes it, separating the
// Gauss-Legendre quadrature error from the interpolation error:
//
//   reference:  the convolution of the exact test function A, integrated to ~1e-12 by a composite quadrature
//               that is independent of Candia's mappings and nodes
//   eps_GL:     Candia's mappings and Gauss-Legendre nodes, but with the exact A at the nodes (no interpolation)
//   eps_total:  what the evolution uses: the precomputed kernel (Grid::buildKernel) acting on A at the grid points
//   eps_interp: the difference of the two Candia results relative to the reference; they use the same nodes,
//               so this is the interpolation error alone, whatever the number of Gauss-Legendre points
//
// Sweep A varies the number of Gauss-Legendre points at the default grid; Sweep B varies the grid density at a
// large number of Gauss-Legendre points (so eps_total is close to eps_interp until the quadrature error is reached).
//
// usage: numerical-precision   (no arguments; writes sweepA-<kernel>.dat and sweepB-<kernel>.dat, see plot.gplt)

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"
using namespace Candia2;

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <format>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace
{
	using test_function = double(*)(double);

	// x times the Les Houches toy gluon and u valence distributions at Q0
	double gluonLike(double x) { return 1.7*std::pow(x, -0.1)*std::pow(1.0-x, 5); }
	double valenceLike(double x) { return 5.1072*std::pow(x, 0.8)*std::pow(1.0-x, 3); }

	struct Case final
	{
		std::string name;
		std::function<std::unique_ptr<Expression>()> make;
		test_function A;
	};

	// one x value in each region of the grid: logarithmic (two), linear and quadratic
	const std::vector<double> xs{1e-4, 1e-2, 0.3, 0.9};
	const std::vector<uint> sweepA_points{5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 500};
	const std::vector<double> sweepB_scales{0.4, 0.5, 0.6, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0};
	constexpr uint sweepB_points = 200;
	const GridFillerArgs default_filler{};

	/** @brief nodes and weights of an n-point Gauss-Legendre rule on [-1, 1] */
	void gaussLegendre(int n, std::vector<double>& X, std::vector<double>& W)
	{
		X.resize(n);
		W.resize(n);
		for (int i=0; i<n; ++i) {
			double z = std::cos(PI*(i+0.75)/(n+0.5)), step{};
			auto legendre = [n](double z, double& p1, double& p2) {
				p1 = 1.0;
				p2 = 0.0;
				for (int j=1; j<=n; ++j) {
					const double p3 = p2;
					p2 = p1;
					p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j;
				}
			};
			double p1{}, p2{};
			do {
				legendre(z, p1, p2);
				step = p1/(n*(z*p1 - p2)/(z*z - 1.0));
				z -= step;
			} while (std::abs(step) > 1e-15);
			legendre(z, p1, p2);
			const double pp = n*(z*p1 - p2)/(z*z - 1.0);
			X[i] = z;
			W[i] = 2.0/((1.0 - z*z)*pp*pp);
		}
	}

	/** @brief the integrand of the convolution at x, as in Grid::mappingFunctionBase, with the exact A */
	double integrand(Expression& E, test_function A, double x, double y, double one_minus_y)
	{
		return E.calcRegular(y)*A(x/y) + E.plus()*(A(x/y) - A(x))/one_minus_y;
	}

	/**
	 *  @brief the reference value of the convolution: 20-point Gauss-Legendre on 4000 panels in each of
	 *  y = x e^u for y in [x, 1/2] (resolves the small-y structure) and y = 1 - (1-y0) e^{-t} for y in [y0, 1]
	 *  (resolves the logarithms at y -> 1), up to 1-y = 1e-14, where y = 1 in double precision
	 */
	double reference(Expression& E, test_function A, double x)
	{
		static std::vector<double> X, W;
		if (X.empty())
			gaussLegendre(20, X, W);
		constexpr int panels = 4000;
		long double sum = 0.0L;

		const double ymid = std::max(x, 0.5);
		if (x < ymid) {
			const double du = std::log(ymid/x)/panels;
			for (int p=0; p<panels; ++p) {
				for (size_t i=0; i<X.size(); ++i) {
					const double y = x*std::exp((p + 0.5*(X[i] + 1.0))*du);
					sum += 0.5*du*W[i]*y*integrand(E, A, x, y, 1.0-y);
				}
			}
		}
		const double dt = std::log((1.0-ymid)/1e-14)/panels;
		for (int p=0; p<panels; ++p) {
			for (size_t i=0; i<X.size(); ++i) {
				const double one_minus_y = (1.0-ymid)*std::exp(-(p + 0.5*(X[i] + 1.0))*dt);
				sum += 0.5*dt*W[i]*one_minus_y*integrand(E, A, x, 1.0-one_minus_y, one_minus_y);
			}
		}
		return static_cast<double>(sum) + A(x)*(E.plus()*std::log1p(-x) + E.delta());
	}

	/**
	 *  @brief Candia's quadrature with the exact A at the nodes; the mappings are those of Grid::setupMappings
	 *  and Grid::getMappings (keep in sync), the nodes and weights are the grid's own
	 */
	double quadratureOnly(Grid const& grid, GridFillerArgs const& filler, Expression& E, test_function A, double x)
	{
		const double p1 = filler.pivot1, p2 = filler.pivot2;
		auto mapped = [&](auto&& mapping) {
			double sum = 0.0;
			for (size_t i=0; i<grid.abscissae().size(); ++i) {
				auto [y, J] = mapping(grid.abscissae()[i]);
				sum += grid.weights()[i]*J*integrand(E, A, x, y, 1.0-y);
			}
			return sum;
		};
		auto quadratic = [](double a) {
			return [a](double z) { return std::pair{1.0-(1.0-a)*(1.0-z)*(1.0-z), 2.0*(1.0-a)*(1.0-z)}; };
		};

		double res = A(x)*(E.plus()*std::log1p(-x) + E.delta());
		if (x < p1) {
			res += mapped([&](double z) { const double a = p1/x; return std::pair{x*std::pow(a, z), x*std::pow(a, z)*std::log(a)}; });
			res += mapped([&](double z) { return std::pair{p1 + (p2-p1)*z, p2-p1}; });
			res += mapped(quadratic(p2));
		} else if (x < p2) {
			res += mapped([&](double z) { return std::pair{x + (p2-x)*z, p2-x}; });
			res += mapped(quadratic(p2));
		} else {
			res += mapped(quadratic(x));
		}
		return res;
	}

	GridFillerArgs scaledFiller(double scale)
	{
		GridFillerArgs filler = default_filler;
		filler.log_size = static_cast<uint>(std::lround(default_filler.log_size*scale));
		filler.lin_size = static_cast<uint>(std::lround(default_filler.lin_size*scale));
		filler.quad_size = static_cast<uint>(std::lround(default_filler.quad_size*scale));
		return filler;
	}

	struct Errors final
	{
		std::vector<double> gl{}, total{}, interp{}; //!< eps_GL, eps_total and eps_interp at each of xs
	};

	Errors errors(
		Expression& E, test_function A, std::vector<double> const& refs, GridFillerArgs const& filler, uint gauss_points, uint& grid_size)
	{
		Grid grid(xs, filler, ConvIntArgs{.num_gauss_points=gauss_points, .num_interp_points=4});
		grid_size = grid.size();
		std::vector<double> a(grid.size());
		for (uint k=0; k<grid.size(); ++k)
			a[k] = A(grid[k]);
		const ConvolutionKernel kernel = grid.buildKernel(E);

		Errors out{};
		for (size_t i=0; i<xs.size(); ++i) {
			const uint k = grid.ntab()[i];
			const double quadrature = quadratureOnly(grid, filler, E, A, grid[k]);
			const double total = kernel.apply(a, k);
			out.gl.push_back(std::abs(quadrature/refs[i] - 1.0));
			out.total.push_back(std::abs(total/refs[i] - 1.0));
			out.interp.push_back(std::abs((total - quadrature)/refs[i]));
		}
		return out;
	}

	/** @brief least-squares slope of log(eps) against log(h), over the points well above the precision floor */
	double fittedSlope(std::vector<double> const& h, std::vector<double> const& eps)
	{
		const double floor = *std::min_element(eps.begin(), eps.end());
		double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
		for (size_t i=0; i<h.size(); ++i) {
			if (eps[i] < 30.0*floor)
				continue;
			const double lx = std::log(h[i]), ly = std::log(eps[i]);
			n += 1; sx += lx; sy += ly; sxx += lx*lx; sxy += lx*ly;
		}
		if (n < 3)
			return std::nan("");
		return (n*sxy - sx*sy)/(n*sxx - sx*sx);
	}
}

int main()
{
	getLogOptions().verbosity = LOG_WARNING;

	auto splitting = [](auto make) {
		return [make]() -> std::unique_ptr<Expression> {
			std::unique_ptr<Expression> E = make();
			if (auto* sf = dynamic_cast<SplittingFunction*>(E.get()))
				sf->update(4, 25.0/3.0, 0.0); // n_f = 4, beta0 for n_f = 4, mu_R = mu_F
			E->preCalc();
			return E;
		};
	};
	const std::vector<Case> cases{
		{"p0qq", splitting([]{ return std::make_unique<P0qq>(); }), valenceLike},
		{"p0gg", splitting([]{ return std::make_unique<P0gg>(); }), gluonLike},
		{"p1qq", splitting([]{ return std::make_unique<P1qq>(); }), valenceLike},
		{"p1gg", splitting([]{ return std::make_unique<P1gg>(); }), gluonLike},
		{"p2gg", splitting([]{ return std::make_unique<P2gg>(); }), gluonLike},
		{"p3gg", splitting([]{ return std::make_unique<P3gg>(); }), gluonLike},
	};

	for (Case const& c : cases) {
		std::unique_ptr<Expression> E = c.make();
		std::vector<double> refs{};
		for (double x : xs)
			refs.push_back(reference(*E, c.A, x));

		// Sweep A: number of Gauss-Legendre points at the default grid
		{
			std::ofstream out(std::format("sweepA-{}.dat", c.name));
			out << "# Sweep A for " << c.name << ": default grid; columns: N_GL, then eps_GL and eps_total at x =";
			for (double x : xs)
				out << ' ' << x;
			out << '\n';
			for (uint n : sweepA_points) {
				uint size{};
				const Errors e = errors(*E, c.A, refs, default_filler, n, size);
				out << n;
				for (size_t i=0; i<xs.size(); ++i)
					out << std::format(" {:.4e} {:.4e}", e.gl[i], e.total[i]);
				out << '\n';
			}
		}

		// Sweep B: grid density at many Gauss-Legendre points
		std::vector<double> h_rel{};
		std::vector<std::vector<double>> eps(xs.size());
		{
			std::ofstream out(std::format("sweepB-{}.dat", c.name));
			out << "# Sweep B for " << c.name << ": N_GL = " << sweepB_points << "; columns: h/h_default (= 1/scale), grid points,"
				<< " h_log (decades per point), h_lin, h_quad (spacing of sqrt((1-x)/(1-pivot2))), then eps_total and eps_interp at x =";
			for (double x : xs)
				out << ' ' << x;
			out << '\n';
			for (double scale : sweepB_scales) {
				const GridFillerArgs filler = scaledFiller(scale);
				uint size{};
				const Errors e = errors(*E, c.A, refs, filler, sweepB_points, size);
				const double decades = std::round(std::log10(filler.pivot1/filler.min));
				const double h_log = 1.0/static_cast<double>(filler.log_size/static_cast<uint>(decades));
				const double h_lin = (filler.pivot2 - filler.pivot1)/filler.lin_size;
				const double h_quad = 1.0/filler.quad_size;
				out << std::format("{:.5f} {} {:.5e} {:.5e} {:.5e}", 1.0/scale, size, h_log, h_lin, h_quad);
				for (size_t i=0; i<xs.size(); ++i) {
					out << std::format(" {:.4e} {:.4e}", e.total[i], e.interp[i]);
					eps[i].push_back(e.interp[i]);
				}
				out << '\n';
				h_rel.push_back(1.0/scale);
			}
		}

		std::printf("%s: Sweep B slope of eps_interp vs h (above the precision floor):", c.name.c_str());
		for (size_t i=0; i<xs.size(); ++i)
			std::printf("  x=%g: %.5f", xs[i], fittedSlope(h_rel, eps[i]));
		std::printf("\n");
	}
}
