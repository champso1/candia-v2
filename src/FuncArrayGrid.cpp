#include "Candia-v2/FuncArrayGrid.hpp"

#include <cassert>
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
		int idx = grid.interpFindIdx(x);
	    base_type proj{_base.begin()+idx, _base.begin()+idx+2*INTERP_POINTS};

		assert(proj.size() == 2*INTERP_POINTS);
		
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		double const* xa = &(grid.points().data()[idx]);
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
					std::println("        m,i=({},{}), ho={}, hp={}, xa[i]={}, xa[i+m]={}, x={}", m, i, ho, hp, xa[i], xa[i+m], x);
					std::println("Grid points: {}", _grid.get().points());
					exit(EXIT_FAILURE);
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


	double ArrayGrid::operator[](uint idx) const
	{
		return _base[idx];
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


	bool FunctionGrid::_split_n3lo_int = false;

	FunctionGrid::FunctionGrid(Grid const& grid, std::unique_ptr<Expression> expr)
		: _grid{std::cref(grid)}, _func{std::move(expr)}
	{
		addFunctionPoints(_grid.get().points());
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
				_func->regular(x),
				_func->plus(x),
				_func->delta(x),
			});
	}


	double FunctionGrid::convolution(ArrayGrid & A, uint k)
	{
		double x = _grid.get().at(k);
		double logx =  std::log(x);

		double res = (_func->plus(1.0)*std::log1p(-x) + _func->delta(1.0)) * A[k];

		if (!_split_n3lo_int)
		{
			for (uint i=0; i<GAUSS_POINTS; i++)
			{
				double y = _grid.get().abscissae(i);
				double w = _grid.get().weights(i);

				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);

				double interp1 = A(b);
				double interp2 = A(a);

				res -= w*logx*a*_func->regular(a)*interp1;
				res -= w*logx*b*(_func->plus(b)*interp2 - _func->plus(1.0)*A[k])/(1.0-b);
			}
		}
		else
		{
			std::vector<double> const& Xi1{_grid.get().abscissae(1)};
			std::vector<double> const& Wi1{_grid.get().weights(1)};
			for (uint i=0; i<GAUSS_POINTS; i++)
			{
				double y = Xi1[i];
				double w = Wi1[i];

				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);

				double interp1 = A(b);
				double interp2 = A(a);

				res -= w*logx*a*_func->regular(a)*interp1;
				res -= w*logx*b*(_func->plus(b)*interp2 - _func->plus(1.0)*A[k])/(1.0-b);
			}

			std::vector<double> const& Xi2{_grid.get().abscissae(2)};
			std::vector<double> const& Wi2{_grid.get().weights(2)};
			for (uint i=0; i<GAUSS_POINTS; i++)
			{
				double y = Xi2[i];
				double w = Wi2[i];

				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);

				double interp1 = A(b);
				double interp2 = A(a);

				res -= w*logx*a*_func->regular(a)*interp1;
				res -= w*logx*b*(_func->plus(b)*interp2 - _func->plus(1.0)*A[k])/(1.0-b);
			}

			std::vector<double> const& Xi3{_grid.get().abscissae(3)};
			std::vector<double> const& Wi3{_grid.get().weights(3)};
			for (uint i=0; i<GAUSS_POINTS; i++)
			{
				double y = Xi3[i];
				double w = Wi3[i];

				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);

				double interp1 = A(b);
				double interp2 = A(a);

				res -= w*logx*a*_func->regular(a)*interp1;
				res -= w*logx*b*(_func->plus(b)*interp2 - _func->plus(1.0)*A[k])/(1.0-b);
			}
		}
			
		return res;
	}


	
}
