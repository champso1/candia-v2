#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <iostream>
#include <limits>
#include <cmath>


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
	}

	void Grid::InitGauLeg()
	{
		const double eps = 1.0e-14; // relative precision

		// abscissae are symmetric:
		uint n = GAUSS_POINTS; // simpler to type
		uint m = (n+1)/2;
		double x2 = _grid_points.at(this->Size()-1);
		double x1 = _grid_points.at(0);
		double xm = 0.5*(x2+x1),
			   xl = 0.5*(x2-x1);

		for (uint i=0; i<m; i++)
		{
			double z = std::cos(PI*(i+0.75)/(n+0.5));

			// default initialize some of these.
			// is a nice check for if they error
			double z1 = std::numeric_limits<double>::max();
			double pp = std::numeric_limits<double>::max();

			double p1, p2, p3;
			double J;
			do
			{
				p1 = 1.0;
				p2 = 0.0;

				for (uint j=0; j<n; j++)
				{
					J = static_cast<double>(j);
					p3 = p2;
					p2 = p1;
					p1 = ((2.0*J + 1.0)*z*p2 - J*p3)/(J+1.0);
				}

				pp = static_cast<double>(n)*(z*p1 - p2)/(z*z - 1.0);
				z1 = z;
				z = z1 - p1/pp;
				
			} while (std::abs(z-z1) > eps);

			if (z1 == std::numeric_limits<double>::max() ||
				pp == std::numeric_limits<double>::max())
			{
				throw("[MATH] error determining gauss-legendre abscissae/weights");
			}

			_Xi[i] = xm - xl*z;
			_Xi[n-1-i] = xm + xl*z;
			_Wi[i] = 2.0*xl/((1.0 - z*z)*pp*pp);
			_Wi[n-1-i] = _Wi[i];
		}
	}

	Grid::Grid(std::vector<double> const& xtab, const uint nx)
		: _grid_points(nx), _Xi(GAUSS_POINTS), _Wi(GAUSS_POINTS)
	{
		InitGrid(xtab, nx);
		InitGauLeg();
	}

	static void __check_idx(const uint idx, const uint size)
	{
		if (idx >= size)
		{
			std::cerr << "[GRID] Index " << idx << "out of range for size " << size << '\n';
			exit(1);
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

}
