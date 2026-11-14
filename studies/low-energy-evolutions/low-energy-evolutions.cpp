#include <string>
#include <fstream>
#include <vector>
#include <limits>
#include <iomanip>

#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"
using namespace Candia2;
using out_type = std::vector<ArrayGrid>;

/*
  0.0, 0.0, 0.0, 0.0, 1.45, 4.75, 172.0, 0.0
  as(q0=1.0) = 0.475
 */

class LesHouchesDistribution2 final : public Distribution
{
public:
	using typename Distribution::value_type;
	using typename Distribution::accessor_type;
private:
	static constexpr std::array<value_type,8> _modified_masses = { 0.0, 0.0, 0.0, 0.0, 1.45, 4.75, 172.0, 0.0 };
public:
	explicit inline LesHouchesDistribution2(double qf)
	{
		setup(1.0, qf);
	}

	inline value_type xuv(value_type x) const
	{
		return 5.1072*std::pow(x, 0.8)*std::pow(1.0-x, 3.0);
	}
	inline value_type xdv(value_type x) const
	{
		return 3.06432*std::pow(x, 0.8)*std::pow(1.0-x, 4.0);
	}

	inline value_type xg (value_type x) const override
	{
		return 1.7*std::pow(x, -0.1)*std::pow(1.0-x, 5.0);
	}
	inline value_type xu (value_type x) const override
	{
		return xuv(x) + xub(x);
	}
	inline value_type xd (value_type x) const override
	{
		return xdv(x) + xdb(x);
	}
	inline value_type xdb(value_type x) const override
	{
		return 0.1939875*std::pow(x, -0.1)*std::pow(1.0-x, 6.0);
	}
	inline value_type xub(value_type x) const override
	{
		return xdb(x)*(1.0-x);
	}
	inline value_type xs (value_type x) const override
	{
		return 0.2*(xub(x) + xdb(x));
	}
	inline value_type xsb(value_type x) const override { return xs(x); }

	void fillSingletCoeffs(
		accessor_type const& accessor,
		std::vector<value_type> const& grid_points) const override
	{
		for (uint k=0; k<grid_points.size(); k++) {
			double x = grid_points[k];
			accessor(0, k) = xg(x);
			accessor(1, k) = xqplus(x);
		}
	}
	void fillNonSingletCoeffs(
		accessor_type const& accessor,
		std::vector<value_type> const& grid_points) const override
	{
		for (uint k=0; k<grid_points.size(); k++) {
			double x = grid_points[k];
			accessor(1, k) = xu(x);  // u
			accessor(2, k) = xd(x);  // d
			accessor(3, k) = xs(x);  // s
			accessor(7, k) = xub(x); // ub
			accessor(8, k) = xdb(x); // db
			accessor(9, k) = xs(x);  // sb ( = s)
		}
	}

	void setup(double q0, double qf) override
	{
		_Q0 = q0;
		_Qf = qf;
		_alpha0 = 0.475;
		_nfi = 3;
		_nff = 4;
		_masses = _modified_masses;
	}
};

using run_output_type = std::pair<std::vector<double>,std::vector<double>>;

static run_output_type run(uint order, double qf);
static void write_output(
	double energy, Grid const& grid,
    run_output_type const& nnlo,
	run_output_type const& n3lo);

int main()
{
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab, {.min=1e-7, .log_size=150}, {.num_gauss_points=100});
	
	for (double qf : {2.0, 3.0}){
		auto nnlo = run(2, qf);
		auto n3lo = run(3, qf);
		write_output(qf, grid, nnlo, n3lo);
	}
}

static run_output_type run(uint order, double Qf)
{
	std::vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
	Grid grid(xtab, {.min=1e-7, .log_size=150}, {.num_gauss_points=100});
	
	std::string datafile_name = std::format("n{}lo-qf{}.dat", order, Qf);
	constexpr uint iterations = 20;
	constexpr uint trunc_idx = 15;
	constexpr double mur2_muf2 = 1.0;

	LesHouchesDistribution2 dist(Qf);
	AlphaS alphas(order, dist.Q0(), Qf, dist.alpha0(), mur2_muf2);
	alphas.setVFNS(dist.masses(), dist.nfi(), dist.nff());
	
	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, dist, mur2_muf2);

	auto F = solver.evolve();
	std::vector<double> c{};
	std::vector<double> g{};
	for (uint k : grid.ntab()) {
		c.emplace_back(F[4][k]);
		g.emplace_back(F[0][k]);
	}
	return std::make_pair(c, g);
}

static void write_output(
    double energy, Grid const& grid,
    run_output_type const& nnlo,
	run_output_type const& n3lo)
{
	std::string datafile_name = std::format("qf{}.dat", energy);
	std::ofstream outfile(datafile_name);
	outfile << std::scientific;
	for (uint i=0; i<grid.xtab().size(); ++i) {
		double x = grid.xtab()[i];
		double nnlo_c = nnlo.first[i];
		double nnlo_g = nnlo.second[i];
		double n3lo_c = n3lo.first[i];
		double n3lo_g = n3lo.second[i];

		auto calc = [](double nnlo_val, double n3lo_val) {
			return (n3lo_val-nnlo_val)/nnlo_val * 100.0;
		};
		
		outfile << std::setprecision(2) << x << ' ';
		outfile << std::setprecision(std::numeric_limits<double>::max_digits10);
		outfile << calc(nnlo_c, n3lo_c) << ' ' << calc(nnlo_g, n3lo_g) << '\n';
	}
}
