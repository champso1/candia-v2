#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <cstdlib>
#include <random>
#include <iostream>
#include <iterator>
#include <limits>
#include <cmath>
#include <memory>
#include <set>


namespace Candia2
{
	Grid::Grid(std::vector<double> const& xtab, const uint nx)
		: _points(nx), _Xi(GAUSS_POINTS), _Wi(GAUSS_POINTS)
	{
		InitGrid(xtab, nx);
		InitGauLeg();
	}
	
	void Grid::InitGrid(std::vector<double> const& xtab, const uint nx)
	{
		std::cerr << "[GRID] Using InitGrid()\n";
		
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
				_points.at(j) = xtab[i]*std::pow(10.0, lstep*(double)(j-ntab[i]));
		}

		_points.at(nx-1) = 1.0;
		_ntab = ntab;
	}

	void Grid::InitGrid2(std::vector<double> const& xtab, const uint nx)
	{
		std::cerr << "[GRID] Using InitGrid2()\n";
		
		auto tanh_sinh_map = [](double t) -> double
		{
			return 0.5 * (1.0 + std::tanh((M_PI / 2.0) * std::sinh(t)));
		};

		// fill the array with the grid points
		// fill with nx-xtab.size(),
		// since we will fill those in later
		static const double T = 3.0;
		uint size = nx - xtab.size();
		std::vector<double> points(size);
		double SIZE = static_cast<double>(size), I{};
		for (uint i=0; i<size; ++i)
		{
			I = static_cast<double>(i);
			double t = -T + (2.0*T)*I / (SIZE-1.0); // linear spacing in t
			points.at(i) = tanh_sinh_map(t);
		}

		// scale them between the min tabulated value and 1 - 1e-15;
		double xmin = xtab.front();
		double xmax = 1.0-1e-5;
		double pmin = points.front();
		double pmax = points.back();
		double xrange = xmax - xmin;
		double prange = pmax - pmin;
		for (double &p : points)
			p = xmin + (xrange/prange)*(p-pmin);

		
		// insert the tabulated points
		std::ranges::sort(points);
		std::set<double> set{points.begin(), points.end()};
		set.insert(xtab.begin(), xtab.end());

		// check to make sure we have enough elements
		// if not just add some random points
		std::random_device rd{};
		std::mt19937 gen{rd()};
		std::uniform_real_distribution<double> dis(points.front(), points.back());
		if (set.size() < nx)
		{
			std::cerr << "[GRID] Grid::InitGrid2(): set size insufficient after inserting xtab.\n";
			while (nx - set.size() > 0)
			{
			    set.insert(dis(gen));
			}
		}
		
		_points = std::vector<double>{set.begin(), set.end()};
		std::ranges::sort(_points);

		// build the ntab array
		for (const double x : xtab)
		{
			auto it = std::ranges::lower_bound(_points, x);
			if (it != _points.end() && std::abs(*it - x) < 1e-14)
				_ntab.emplace_back(std::distance(_points.begin(), it));
		}

		// set the absolute last value to exactly 1
		_points.back() = 1.0;

		auto def_precision = std::cerr.precision();
		std::cerr.precision(std::numeric_limits<double>::max_digits10);
		std::cerr << "[GRID] Grid::InitGrid2(): Points array is: [\n";
		for (const double n : _points)
			std::cerr << "    " << n << '\n';
		std::cout << "]\n";
		std::cerr.precision(def_precision);
		
		std::cerr << "[GRID] Grid::InitGrid2(): Ntab array is: [";
		for (const uint n : _ntab)
			std::cerr << ' ' << n << ' ';
		std::cout << "]\n";
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
				std::cerr << "[ERROR] Grid::InitGauleg(): error determining gauss-legendre abscissae/weights\n";
				exit(EXIT_FAILURE);
			}

			_Xi[i-1] = xm - xl*z;
			_Xi[n-i] = xm + xl*z;
			_Wi[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp);
			_Wi[n-i] = _Wi[i-1];
		}
	}

    uint Grid::InterpFindIdx(double x) const
	{
		int k;
		for (k=0; x>=_points.at(k); k++)
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

	

	double Grid::Interpolate(std::vector<double> const& yy, const double x) const
	{
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = static_cast<int>(this->InterpFindIdx(x));

		double const* xa = &(_points.data()[k]);
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

				den = ho-hp;
				if (std::abs(ho-hp) < 1e-15)
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
		double x = _points.at(k);
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




	double FunctionGrid::Interpolate(FunctionPart part, const double x) const
	{
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = static_cast<int>(_grid.InterpFindIdx(x));

		double const* xa = &(_grid.Points().data()[k]);
		double const* ya = &(Y(part).data()[k]);

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

	double FunctionGrid::Convolution(std::vector<double> A, uint k) const
	{
		double x = _grid[k];
		double logx =  std::log(x);

		double eplus = Y(PLUS, 1.0);
		double edelta = Y(DELTA, 1.0);
		double res = (eplus*std::log1p(-x) + edelta) * A.at(k);

		for (uint i=0; i<GAUSS_POINTS; i++)
		{
		    double y = _grid.Abscissae(i);
			double w = _grid.Weights(i);

			double a = std::pow(x, 1.0-y);
			double b = std::pow(x, y);

			double interp1 = _grid.Interpolate(A, b);
			double interp2 = _grid.Interpolate(A, a);

			res -= w*logx*a*Interpolate(REGULAR, a)*interp1;
			res -= w*logx*b*(Interpolate(PLUS, b)*interp2 - eplus*A.at(k))/(1.0-b);
		}
		
		return res;
	}

}
