#include "Candia-v2/FuncArrGrid.hpp"

#include <cassert>
#include <functional>
#include <cmath>
#include <print>

namespace Candia2
{
	ArrayGrid::size_type ArrayGrid::size() const noexcept
	{
		return _cache.size();
	}

	void ArrayGrid::zero() noexcept
	{
		_cache.clear();
		for (double& _x : _base)
			_x = 0.0;
	}

	double ArrayGrid::interpolate(double x)
	{
		Grid const& grid = _grid.get();
		int idx = grid.InterpFindIdx(x);
	    base_type proj{_base.begin()+idx, _base.begin()+idx+2*INTERP_POINTS};

		assert(proj.size() == 2*INTERP_POINTS);
		
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		double const* xa = &(grid.Points().data()[idx]);
		double const* ya = proj.data();
		
	    base_type c(n, 0.0);
	    base_type d(n, 0.0);

		dif = std::abs(x - xa[0]);

		for (int i=0; i<n; i++)
		{
			if ((dift = std::abs(x - xa[i])) < dif)
			{
				ns = i;
				dif = dift;
			}
			c[i] = ya[i];
			d[i] = ya[i];
		}

		y = ya[ns--];

		for (int m=1; m<n; m++)
		{
			for (int i=0; i<n-m; i++)
			{
				ho = xa[i] - x;
				hp = xa[i+m] - x;
				w = c[i+1] - d[i];

				den = ho-hp;
				if (std::abs(ho-hp) < 1e-15)
				{
					std::println("[AGRID] Interpolate(): found a denominator equal to 0.0.");
					exit(1);
				}

				den = w/den;
				d[i] = hp*den;
				c[i] = ho*den;
			}

			y += (2*ns < (n-1-m) ? c[ns+1] : d[ns--]);
		}
			
		return y;
	}


	void ArrayGrid::addPoint(double x)
	{
		double val = interpolate(x);
		_cache[x] = val;
	}


	double& ArrayGrid::operator[](uint idx)
	{
		return _base[idx];
	}

	double& ArrayGrid::operator()(double x)
	{
		if (_cache.find(x) != _cache.end())
		{
			return _cache.at(x);
		}

		addPoint(x);
		return _cache.at(x);
	}



	FunctionGrid::FunctionGrid(Grid const& grid, std::unique_ptr<Expression> expr)
		: _grid{std::cref(grid)}, _func{std::move(expr)}
	{
		addFunctionPoints(_grid.get().Points());
	}

	double FunctionGrid::operator()(double x, uint function_part)
	{
		if (_cache.find(x) != _cache.end())
		{
			return _cache.at(x).at(function_part);
		}

		addFunctionPoint(x);
		return _cache.at(x).at(function_part);
	}

	void FunctionGrid::addFunctionPoints(std::vector<double> const& X)
	{
		for (double x : X)
			addFunctionPoint(x);
	}

	void FunctionGrid::addFunctionPoint(double x)
	{
		_cache.emplace(x, std::array<double, 3>{
				_func->Regular(x),
				_func->Plus(x),
				_func->Delta(x),
			});
	}


	double FunctionGrid::convolution(ArrayGrid & A, uint k)
	{
		double x = _grid.get().At(k);
		double logx =  std::log(x);

		double res = (_func->Plus(1.0)*std::log1p(-x) + _func->Delta(1.0)) * A[k];

		for (uint i=0; i<GAUSS_POINTS; i++)
		{
			double y = _grid.get().Abscissae(i);
			double w = _grid.get().Weights(i);

			double a = std::pow(x, 1.0-y);
			double b = std::pow(x, y);

			double interp1 = A(b);
			double interp2 = A(a);

			res -= w*logx*a*_func->Regular(a)*interp1;
			res -= w*logx*b*(_func->Plus(b)*interp2 - _func->Plus(1.0)*A[k])/(1.0-b);
		}
			
		return res;
	}


	
}
