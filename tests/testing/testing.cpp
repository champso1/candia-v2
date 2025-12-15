#include "Candia-v2/Grid.hpp"
#include "Candia-v2/OperatorMatrixElements.hpp"
using namespace Candia2;

#include <map>
#include <memory>
#include <concepts>
using namespace std;

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

int main(int argc, char *argv[])
{
	(void)argc;
	(void)argv;

	const uint num_grid_points = 100;
	vector<double> xtab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, num_grid_points);

	createExpression<OpMatElemN3LO>("A3psqq", grid, ome::AqqQPS);
	OpMatElemN3LO::update(0, 3);
	for (auto& [name, expr] : _expressions)
		expr->fill(grid.points(), grid.abscissae());
	Expression& expr = getExpression("A3psqq");

	auto _a3psqq = ome::AqqQPS;
	auto a3psqq = _a3psqq.get_regular().value();

	double test_point = std::pow(grid[10], 1.0-grid.abscissae(10));

	double val1 = a3psqq[3](0, 3, test_point);
	double val2 = expr.regular(test_point);

	std::println("val1={}", val1);
	std::println("val2={}", val2);
}
