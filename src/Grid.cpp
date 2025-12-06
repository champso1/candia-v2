#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <set>
#include <print>


namespace Candia2
{
	Grid::Grid(std::vector<double> const& xtab, uint nx, int grid_fill_type)
		: _points(nx), _ntab{},
		  _Xi(GAUSS_POINTS), _Wi(GAUSS_POINTS),
		  _Xi1(GAUSS_POINTS), _Wi1(GAUSS_POINTS),
		  _Xi2(GAUSS_POINTS), _Wi2(GAUSS_POINTS),
		  _Xi3(GAUSS_POINTS), _Wi3(GAUSS_POINTS)
	{
		switch (grid_fill_type)
		{
			case 1: initGrid(xtab, nx); break;
			case 2: initGrid2(xtab, nx); break;
			case 3: initGrid3(xtab, nx); break;
			default:
			{
				std::println(stderr, "[DGLAP2] Warning: Invalid grid fill type. Found {}, expected 1, 2, or 3.", grid_fill_type);
				std::println(stderr, "         Will use default (1).");
				initGrid(xtab, nx);
			}
		}
	    
		initGauLeg(0.0, 1.0, _Xi, _Wi);

		initGauLeg(0, 1e-1, _Xi1, _Wi1);
		initGauLeg(1e-1, 0.7, _Xi2, _Wi2);
		initGauLeg(0.7, 1.0, _Xi3, _Wi3);
	}

	void Grid::initGrid(std::vector<double> const& xtab, const uint nx)
	{
		std::println("[GRID] Using InitGrid()");
		
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

		return;
	}


	void Grid::initGrid2(std::vector<double> const& xtab, uint nx)
	{
		std::println("[GRID] Using InitGrid2()");
		
		std::vector<double> points(nx-xtab.size()+1);
		const double xmin = xtab.front();
		const double xmax = xtab.back();
		const double ymin = std::log10(xmin);
		const double ymax = std::log10(xmax);
		for (uint i=0; i<nx-xtab.size()+1; ++i)
		{
			const double u = static_cast<double>(i)/(static_cast<double>(nx-xtab.size()+1) - 1.0);
			// double _y = ymin + (ymax-ymin)*0.5*(1.0 - std::cos(PI*u));
			double _y = ymin + (ymax-ymin)*u;
			points[i] = std::pow(10.0, _y);
		}

		// insert the tabulated points
		// make it a set to avoid duplicate values
		// then replace the original points array with the new one
		std::ranges::sort(points);
		std::set<double> set{points.begin(), points.end()};
		set.insert(xtab.begin(), xtab.end());
		points = std::vector<double>(set.begin(), set.end());
		std::ranges::sort(points);

		_points = points;

		// build the ntab array
		_ntab = std::vector<int>{};
		for (const double x : xtab)
		{
			auto it = std::ranges::find(_points, x);
			if (it == _points.end())
			{
				std::println(stderr, "[GRID] Somehow found a tabulated value ({}) that is not in the ntab array.", x);
				exit(EXIT_FAILURE);
			}
			_ntab.emplace_back(std::distance(_points.begin(), it));
		}
	}

	void Grid::initGrid3(std::vector<double> const& xtab, uint nx)
	{
		std::println("[GRID] Using InitGrid3()");
		
		std::vector<double> points{};

		std::vector<double> log_tab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.7, 1.0};
		std::vector<double> log_xtab{log_tab};
		std::ranges::transform(log_xtab, log_xtab.begin(), [](double x) -> double{ return std::log10(x); });
		int num_grid_points_per_bin = nx / xtab.size();

		for (uint i=0; i<log_xtab.size()-1; ++i)
		{
			double logmin = log_xtab[i];
			double logmax = log_xtab[i+1];
			double dlog = (logmax-logmin)/static_cast<double>(num_grid_points_per_bin);

			int num = 0;
			for (double _l=logmin; _l<logmax && num<num_grid_points_per_bin; _l+=dlog, ++num)
				points.emplace_back(std::pow(10, _l));
		}
		
		// insert the tabulated points
		// make it a set to avoid duplicate values
		// then replace the original points array with the new one
		std::ranges::sort(points);
		std::set<double> set{points.begin(), points.end()};
		set.insert_range(xtab);
		points = std::vector<double>(set.begin(), set.end());
		
		_points = points;

		// build the ntab array
		_ntab = std::vector<int>{};
		for (const double x : xtab)
		{
			auto it = std::ranges::lower_bound(_points, x);
			if (it != _points.end() && std::abs(*it - x) < 1e-14)
				_ntab.emplace_back(std::distance(_points.begin(), it));
		}
	}

	void Grid::initGauLeg(double x1, double x2, std::vector<double> & Xi, std::vector<double> & Wi)
	{
		const double eps = 3.0e-11; // relative precision

		// abscissae are symmetric:
		uint n = GAUSS_POINTS; // simpler to type
		double N = static_cast<double>(n);
		uint m = (n+1)/2;
		// double x2 = 1.0;
		// double x1 = 0.0;
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
				std::println("[GRID: ERROR] InitGauleg(): failed to determine gauss-legendre abscissae/weights");
				exit(EXIT_FAILURE);
			}

			Xi[i-1] = xm - xl*z;
			Xi[n-i] = xm + xl*z;
			Wi[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp);
			Wi[n-i] = Wi[i-1];
		}
	}

    uint Grid::interpFindIdx(double x) const
	{
		int k;
		for (k=0; x>=_points.at(k); k++)
		{
			if (k >= static_cast<int>(size()-1))
				break;
		}

		k-=INTERP_POINTS;

		if (k<0) k=0;
		if (k>static_cast<int>(size()-2*INTERP_POINTS))
			k=size()-2*INTERP_POINTS;

		return static_cast<uint>(k);
	}

	

	double Grid::interpolate(std::vector<double> const& yy, const double x) const
	{
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = static_cast<int>(this->interpFindIdx(x));

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
					std::println("[GRID: ERROR] Interpolate(): found a denominator equal to 0.0.");
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



	double Grid::convolution(std::vector<double> const& A,
		std::shared_ptr<Expression> E, uint k)
	{
		double x = _points.at(k);
		double logx =  std::log(x);
		double res = (E->plus(1.0)*std::log1p(-x) + E->delta(1.0)) * A.at(k);
		
		for (uint i=0; i<GAUSS_POINTS; i++)
		{
			double y = _Xi[i];
			double w = _Wi[i];

			double a = std::pow(x, 1.0-y);
			double b = std::pow(x, y);

			double interp1 = interpolate(A, b);
			double interp2 = interpolate(A, a);

			res -= w*logx*a*E->regular(a)*interp1;
			res -= w*logx*b*(E->plus(b)*interp2 - E->plus(1.0)*A.at(k))/(1.0-b);
		}
		
		return res;
	}


}
