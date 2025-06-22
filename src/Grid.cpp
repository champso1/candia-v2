#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/SplittingFn.hpp"

#include <iostream>
#include <cmath>
#include <span>

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
				_grid_points[j] = xtab[i]*std::pow(10.0, lstep*(double)(j-ntab[i]));
		}

		_grid_points[nx-1] = 1.0;
	}

	void Grid::InitGauLeg()
	{
		const double eps = 1.0e-14; // relative precision

		// abscissae are symmetric:
		uint n = GAUSS_POINTS; // simpler to type
		uint m = (n+1)/2;
		double xf = this->At(this->Size()-1);
		double x0 = this->At(0);
		double xl = 0.5*(xf-x0),
			   xu = 0.5*(xf+x0);

		for (uint i=0; i<m; i++)
		{
			double z = std::cos(M_PI*(i+0.75)/(n+0.5));

			// default initialize some of these.
			// is a nice check for if they error
			double z1 = std::numeric_limits<double>::max();
			double pp = std::numeric_limits<double>::max();

			double p1, p2, p3;
			do
			{
				p1 = 1.0;
				p2 = 0.0;

				for (uint j=0; j<n; j++)
				{
					p3 = p2;
					p2 = p1;
					p1 = ((2.0*j + 1.0)*z*p2 - j*p3)/(j+1);
				}

				pp = n*(z*p1 - p2)/(z*z - 1.0);
				z1 = z;
				z = z1 - p1/pp;
				
			} while (std::abs(z-z1) > eps);

			if (z1 == std::numeric_limits<double>::max() ||
				pp == std::numeric_limits<double>::max())
			{
				throw("[MATH] error determining gauss-legendre abscissae/weights");
			}

			_Xi[i] = xu - xl*z;
			_Xi[n-1-i] = xu + xl*z;
			
			_Wi[i] = 2.0*xl/((1.0 - z*z)*pp*pp);
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
			throw("Index " + std::to_string(idx) + " out of range; must be less than " + std::to_string(size));
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

    uint Grid::InterpFindIdx(std::vector<double> const& vec, double x) const
	{
		if (vec.size() != this->Size())
		{
			std::cerr << "[GRID] InterpFindIdx(): f(x) size (" << vec.size()
					  << ") and x size (" << vec.size()
					  <<  ") do not match.\n";
			exit(1);
		}

		//	std::cerr << "[GRID] InterpFindIdx(): Finding a good index...\n";

		uint n = vec.size(); // shorthand
		uint mm = INTERP_POINTS; // also shorthand
		uint ju = 0, jm = 0, jl = 0;
		bool ascend = (this->At(n-1) >= this->At(0));
		ju = n - 1;

		// std::cerr << "[GRID] InterpFindIdx(): Setup, before iterations\n";

		while (ju - jl > 1)
		{
			jm = (ju+jl) >> 1;
			if ((x >= this->At(jm)) == ascend)
				jl = jm;
			else
				ju = jm;
		}

		// std::cerr << "[GRID] InterpFindIdx(): Done.\n";

		return MAX(0, MIN(n-mm, jl-((mm-2) >> 1)));
	}

	double Grid::Interpolate(std::vector<double> const& y, const double x) const
	{
		// std::cerr << "[GRID] Interpolate(): Interpolating...\n";
		
		uint k = this->InterpFindIdx(y, x);

		// std::cerr << "[GRID] Interpolate(): found a good k value: " << k << '\n';

		uint ns = 0;
		double yy, den, dif, dift, ho, hp, w, dy;

		// these are like std::string_view, but for vectors
		std::span<const double> xa(_grid_points.begin()+k, INTERP_POINTS);
		std::span<const double> ya(y.begin()+k, INTERP_POINTS);

		dif = std::abs(x - xa[0]);
		
		std::vector<double> c(INTERP_POINTS);
		std::vector<double> d(INTERP_POINTS);

		// std::cerr << "[GRID] Interpolate(): Grabbed spans, other stuff before iterations.\n";

		for (uint i=0; i<INTERP_POINTS; i++)
		{
			if ((dift = std::abs(x - xa[i])) < dif)
			{
				ns = i;
				dif = dift;
			}
			c[i] = ya[i];
			d[i] = ya[i];
		}

		yy = ya[ns--];

	    // std::cerr << "[GRID] Interpolate(): Completed first iterations(s).\n";

		for (uint m=1; m<INTERP_POINTS; m++)
		{
		    for (uint i=0; i<INTERP_POINTS-m; i++)
			{
				ho = xa[i] - x;
				hp = xa[i+m] - x;
				w = c[i+1] - d[i];

				if ((den = ho-hp) == 0.0)
				{
					std::cerr << "[GRID] Interp(): found a denominator equal to 0.0.\n";
					exit(1);
				}

				den = w/den;
				d[i] = hp*den;
				c[i] = ho*den;
			}
			yy += (dy = (2*(ns+1) < (INTERP_POINTS-m) ? c[ns+1] : d[ns--]));
		}

		// std::cerr << "[GRID] Interpolate(): Completed second iterations(s). Now done. Found: " << yy << '\n';

		return yy;
	}



	double Grid::Convolution(std::vector<double> const& A,
							 std::shared_ptr<SplittingFunction> P,
							 uint k)
	{
		// std::cerr << "[GRID] Convolution(): Before everything\n";
		
		double x = this->At(k);
		double fac1 = P->Delta(1.0)*A[k];
		double fac2 = P->Plus(1.0)*A[k]*std::log1p(x);
		double fac3 = 0.0, fac4 = 0.0;

	    // std::cerr << "[GRID] Convolution(): After initial factors, before actual convolution\n";
		
		for (uint i=0; i<GAUSS_POINTS; i++)
		{
			// std::cerr << "[GRID] Convolution(): Inside convolution, grabbing weights\n";
			double z = _Xi[i];
			double w = _Wi[i];

			// std::cerr << "[GRID] Convolution(): Inside convolution, doing variable redefinitions\n";
			
			double a = std::pow(x, 1.0-z);
			double b = std::pow(x, z);

			// std::cerr << "[GRID] Convolution(): Inside convolution, doing the weighting "
			//		  << "interpolating, and grabbing splitting function values\n";
			
			fac2 *= w * a*P->Regular(a)*Interpolate(A, b);
			fac3 *= w * b*(P->Plus(a)*Interpolate(A, a) - P->Plus(1.0)*A[k])/(1.0-b);

		    // std::cerr << "[GRID] Convolution(): Inside convolution, after iteration.\n";
		}

		// std::cerr << "[GRID] Convolution(): After actual convolution, before final factors.\n";
		
		double logx = std::log(this->At(k));
		fac2 *= logx;
		fac3 *= logx;

		// std::cerr << "[GRID] Convolution(): After everything\n";
		
		return (fac1 + fac2 + fac3 + fac4);
	}

}
