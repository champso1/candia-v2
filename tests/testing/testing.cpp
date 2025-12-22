#include <map>
#include <memory>
#include <concepts>
#include <print>
using namespace std;

#include "Candia-v2/Candia.hpp"
using namespace Candia2;
using arr_type = vector<ArrayGrid>;

#include <gsl/gsl_integration.h>

std::map<std::string_view, std::unique_ptr<Expression>> _expressions{};

template <typename TExpr, typename... TExprArgs>
requires (std::derived_from<TExpr, Expression>)
void createExpression(std::string_view name, Grid const& grid, TExprArgs&&... args)
{
	std::unique_ptr<Expression> ptr = std::make_unique<TExpr>(std::forward<TExprArgs>(args)...);
	_expressions.emplace(std::make_pair(name, std::move(ptr)));
}
Expression& getExpression(std::string_view name)
{
	return *_expressions[name];
}

static double func(double x, void* params)
{
	LesHouchesDistribution dist = *reinterpret_cast<LesHouchesDistribution*>(params);
	double sigma = dist.xuv(x)
		+ 2.0*dist.xub(x)
		+ dist.xdv(x)
		+ 2.0*dist.xdb(x)
		+ 2.0*dist.xs(x);
	double g = dist.xg(x);
	return sigma + g;
}

static double funcu(double x, void* params)
{
	LesHouchesDistribution dist = *reinterpret_cast<LesHouchesDistribution*>(params);
	double u = dist.xuv(x) + dist.xub(x);
	double ub = dist.xub(x);
	return (u - ub)/x;
}

static double funcd(double x, void* params)
{
	LesHouchesDistribution dist = *reinterpret_cast<LesHouchesDistribution*>(params);
	double d = dist.xdv(x) + dist.xdb(x);
	double db = dist.xdb(x);
	return (d - db)/x;
}


static void test_sum_rules()
{
	LesHouchesDistribution _dist{};
    void * params = reinterpret_cast<void*>(&_dist);
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	{
		gsl_function F;
		F.function = &func;
		F.params = params;
		double res{}, error{};
		gsl_integration_qags(&F, 0.0, 1.0, 0.0, 1e-5, 1000, w, &res, &error);
		println("Testing momentum sum rule: {} +- {}", res, error);
	}
	{
		gsl_function F;
		F.function = &funcu;
		F.params = params;
		double res{}, error{};
		gsl_integration_qags(&F, 0.0, 1.0, 0.0, 1e-5, 1000, w, &res, &error);
		println("Testing valence sum rule (u): {} +- {}", res, error);
	}
	{
		gsl_function F;
		F.function = &funcd;
		F.params = params;
		double res{}, error{};
		gsl_integration_qags(&F, 0.0, 1.0, 0.0, 1e-5, 1000, w, &res, &error);
		println("Testing valence sum rule (d): {} +- {}", res, error);
	}
	gsl_integration_workspace_free(w);

	
}

struct Params
{
	Grid& grid;
	arr_type& arr;
};

static double func_final(double x, void* params_)
{
	Params params = *reinterpret_cast<Params*>(params_);
	double res = 0.0;
	for (uint j=1; j<=5; ++j)
		res += params.grid.interpolate(params.arr[j], x) + params.grid.interpolate(params.arr[j+6], x);
	return res + params.grid.interpolate(params.arr[0], x);
}

static double funcu_final(double x, void* params_)
{
    Params params = *reinterpret_cast<Params*>(params_);
	double u = params.grid.interpolate(params.arr[1], x);
	double ub = params.grid.interpolate(params.arr[1+6], x);
	return (u - ub)/x;
}

static double funcd_final(double x, void* params_)
{
	Params params = *reinterpret_cast<Params*>(params_);
	double d = params.grid.interpolate(params.arr[2], x);
	double db = params.grid.interpolate(params.arr[2+6], x);
	return (d - db)/x;
}

static void test_sum_rules_final(Grid& grid, arr_type& dists)
{
	Params params_{grid, dists};
    void * params = reinterpret_cast<void*>(&params_);
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

	{
		gsl_function F;
		F.function = &func_final;
		F.params = params;
		double res{}, error{};
		gsl_integration_qags(&F, 0.0, 1.0, 0.0, 1e-5, 1000, w, &res, &error);
		println("Testing final momentum sum rule: {} +- {}", res, error);
	}
	{
		gsl_function F;
		F.function = &funcu_final;
		F.params = params;
		double res{}, error{};
		gsl_integration_qags(&F, 0.0, 1.0, 0.0, 1e-5, 1000, w, &res, &error);
		println("Testing final valence sum rule (u): {} +- {}", res, error);
	}
	{
		gsl_function F;
		F.function = &funcd_final;
		F.params = params;
		double res{}, error{};
		gsl_integration_qags(&F, 0.0, 1.0, 0.0, 1e-5, 1000, w, &res, &error);
		println("Testing final valence sum rule (d): {} +- {}", res, error);
	}
	gsl_integration_workspace_free(w);

	
}


int main(int argc, char *argv[])
{
	(void)argc;
	(void)argv;

	const uint order = 3;
	const uint num_grid_points = 800;
	const uint iterations = 10;
	const uint trunc_idx = 15;
	const double kr = 1.0;
	const double Qf = 100.0;
	
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, num_grid_points, 3);

	std::unique_ptr<LesHouchesDistribution> dist = std::make_unique<LesHouchesDistribution>();
	AlphaS alphas(order, dist->Q0(), Qf, dist->alpha0(), kr);
	alphas.setVFNS(dist->masses(), dist->nfi());
	// alphas.setFFNS(4);

	DGLAPSolver solver(order, grid, alphas, Qf, iterations, trunc_idx, *dist, kr, true);
	arr_type arr = solver.evolve();
	test_sum_rules_final(grid, arr);
}
