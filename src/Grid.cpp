#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_math.h"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cmath>
#include <format>
#include <memory>


namespace Candia2
{
	void Grid::InitGrid(std::vector<double> const& xtab, const uint nx)
	{
		const uint xtab_len = xtab.size();
		std::vector<double> Ntab(xtab_len);
	    std::vector<int> ntab(xtab_len);

		double temp = -std::log10(xtab[0]);

		for (uint i=1; i<xtab_len; i++)
			Ntab[i] = (double)(nx-1)*std::log10(xtab[i]/xtab[i-1])/temp;

		ntab[0] = nx-1;

		for (uint i=1; i<xtab_len; i++)
		{
			ntab[i] =  (int)Ntab[i];
			Ntab[i] -= (double)ntab[i];
			ntab[0] -= ntab[i];
		}

		for (uint i=1; i<xtab_len; i++)
		{
			if (ntab[i] == 0)
			{
				ntab[i] =  1;
				Ntab[i] -= 1.0;
				ntab[0] -= 1;
			}
		}

		uint n;
		for ( ; ntab[0]<0; ntab[0]++)
		{
			n=0;

			for (uint i=1; i<xtab_len; i++)
			{
				if (ntab[i] != 1)
				{
					if ((n == 0) || (Ntab[i] <= temp))
					{
						n = i;
						temp = Ntab[i];
					}
				}
			}

			ntab[n]--;
			Ntab[n] += 1.0;
		}

		for ( ; ntab[0]>0; ntab[0]--)
		{
			n=0;

			for (uint i=1; i<xtab_len; i++)
			{
				if ((n == 0) || (Ntab[i] > temp))
				{
					n = i;
					temp = Ntab[i];
				}
			}

			ntab[n]++;
			Ntab[n] -= 1.0;
		}

		for (uint i=1; i<xtab_len; i++)
			ntab[i] += ntab[i-1];

		double lstep;
		for (uint i=0; i<xtab_len-1; i++)
		{
			lstep=std::log10(xtab[i+1]/xtab[i])/(double)(ntab[i+1]-ntab[i]);

			for (int j=ntab[i]; j<ntab[i+1]; j++)
				this->At(j) = xtab[i]*std::pow(10.0, lstep*(double)(j-ntab[i]));
		}

		_grid_points.at(nx-1) = 1.0;
		_ntab = ntab;
	}

	void Grid::InitGauLeg()
	{
		const double eps = 3.0e-11; // relative precision

		// abscissae are symmetric:
		uint n = GAUSS_POINTS; // simpler to type
		double N = static_cast<double>(n);
		uint m = (n+1)/2;
		double x2 = 1.0;
		double x1 = 0.0;
		double xm = 0.5*(x2+x1);
		double xl = 0.5*(x2-x1);

		for (uint i=1; i<=m; i++)
		{
			double I = static_cast<double>(i);
			double z = std::cos(PI*(I-0.25)/(N+0.5));

			// default initialize some of these.
			// easy to spot if error occurs
			double z1 = std::numeric_limits<double>::max();
			double pp = std::numeric_limits<double>::max();

			double p1, p2, p3;
			double J;
			do
			{
				p1 = 1.0;
				p2 = 0.0;

				for (uint j=1; j<=n; j++)
				{
					J = static_cast<double>(j);
					p3 = p2;
					p2 = p1;
					p1 = ((2.0*J - 1.0)*z*p2 - (J-1.0)*p3)/J;
				}

				pp = N*(z*p1 - p2)/(z*z - 1.0);
				z1 = z;
				z = z1 - p1/pp;
			} while (std::abs(z-z1) > eps);

			if (z1 == std::numeric_limits<double>::max() ||
				pp == std::numeric_limits<double>::max())
			{
				std::cerr << "Grid::InitGauleg(): error determining gauss-legendre abscissae/weights\n";
				exit(EXIT_FAILURE);
			}

			_Xi[i-1] = xm - xl*z;
			_Xi[n-i] = xm + xl*z;
			_Wi[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp);
			_Wi[n-i] = _Wi[i-1];
		}
	}

	Grid::Grid(std::vector<double> const& xtab, const uint nx)
		: _grid_points(nx), _Xi(GAUSS_POINTS), _Wi(GAUSS_POINTS),
		  _ws{gsl_integration_workspace_alloc(_N),
		    [](gsl_integration_workspace *p){ gsl_integration_workspace_free(p); }}
	{
		InitGrid(xtab, nx);
		InitGauLeg();
	}

	static void __check_idx(const uint idx, const uint size)
	{
		if (idx >= size)
		{
			std::cerr << std::format("[GRID] Index {} out of range for size {}", idx, size) << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	double const& Grid::operator[](const uint idx) const
	{
		__check_idx(idx, _grid_points.size());
		return _grid_points.at(idx);
	}

	double& Grid::operator[](const uint idx)
	{
		__check_idx(idx, _grid_points.size());
		return _grid_points.at(idx);
	}


    uint Grid::InterpFindIdx(double x) const
	{
		int k;
		for (k=0; x>=_grid_points.at(k); k++)
		{
			if (k >= static_cast<int>(Size()-1))
				break;
		}

		k-=INTERP_POINTS;

		if (k<0) k=0;
		if (k>static_cast<int>(Size()-2*INTERP_POINTS))
			k=Size()-2*INTERP_POINTS;

		return static_cast<uint>(k);
	}

	double Grid::Interpolate(std::vector<double> const& yy, const double x, bool debug) const
	{
		UNUSED(debug);

		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = static_cast<int>(this->InterpFindIdx(x));

		double const* xa = &(_grid_points.data()[k]);
		double const* ya = &(yy.data()[k]);
		std::vector<double> c(n, 0.0);
		std::vector<double> d(n, 0.0);

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

				if ((den = ho-hp) == 0.0)
				{
					std::cerr << "[GRID] Interpolate(): found a denominator equal to 0.0.\n";
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



	double Grid::Convolution(std::vector<double> const& A,
							 std::shared_ptr<Expression> E,
							 uint k)
	{
		double x = _grid_points.at(k);
		double logx =  std::log(x);

		double res = (E->Plus(1.0)*std::log1p(-x) + E->Delta(1.0)) * A.at(k);

		for (uint i=0; i<GAUSS_POINTS; i++)
		{
		    double y = _Xi[i];
			double w = _Wi[i];

			double a = std::pow(x, 1.0-y);
			double b = std::pow(x, y);

			double interp1 = Interpolate(A, b);
			double interp2 = Interpolate(A, a);

			res -= w*logx*a*E->Regular(a)*interp1;
			res -= w*logx*b*(E->Plus(b)*interp2 - E->Plus(1.0)*A.at(k))/(1.0-b);
		}
		
		return res;
	}


	double Grid::ConvolutionGSL(
		std::vector<double> const& A,
		std::shared_ptr<Expression> P,
		uint k)
	{
		double x = _grid_points.at(k);
		// std::cout << x << '\n';
		
		GSLConvObj obj{this, k, x, A, P};
		
		auto func_reg = [](double z, void *params) -> double {
			GSLConvObj *obj = reinterpret_cast<GSLConvObj*>(params);
			Grid const* grid = obj->grid;
			const double x = obj->x;
			std::vector<double> const& A = obj->A;
			const std::shared_ptr<Expression> P = obj->P;

			const double a = std::pow(x, 1.0-z);
			const double b = std::pow(x, z);
			return a*P->Regular(a)*grid->Interpolate(A, b);
		};

		auto func_sing = [](double z, void *params) -> double {
			GSLConvObj *obj = reinterpret_cast<GSLConvObj*>(params);
			Grid const* grid = obj->grid;
			const uint k = obj->k;
			const double x = obj->x;
			std::vector<double> const& A = obj->A;
			const std::shared_ptr<Expression> P = obj->P;

			const double a = std::pow(x, 1.0-z);
			const double b = std::pow(x, z);
			return b*(P->Plus(b)*grid->Interpolate(A, a) - P->Plus(1.0)*A.at(k))/(1.0-b);
		};

		gsl_function F;
		F.function = func_reg;
		F.params = reinterpret_cast<void*>(&obj);

		double res_reg{}, res_sing{};
		double abserr{};
		gsl_integration_workspace *ws = _ws.get();

		gsl_integration_qag(&F, 0.0, 1.0, 1e-4, 0.0, _N, GSL_INTEG_GAUSS61, ws, &res_reg, &abserr);

		F.function = func_sing;
		gsl_integration_qag(&F, 0.0, 1.0, 1e-4, 0.0, _N, GSL_INTEG_GAUSS61, ws, &res_sing, &abserr);

		double res = P->Delta(1.0)*A.at(k);
		res += -std::log(x)*res_reg;
		res += -std::log(x)*res_sing + P->Plus(1.0)*A.at(k)*std::log1p(-x);
		return res;
	}
}
