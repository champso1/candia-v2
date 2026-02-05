#include "Candia-v2/Grid.hpp"
#include "Candia-v2/Common.hpp"
#include "Candia-v2/ArrayGrid.hpp"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <set>

namespace Candia2
{
	Grid::Grid(std::vector<double> const& xtab, uint nx, uint gauss_points, int grid_fill_type)
		: _points(nx), _ntab{}, _xtab{xtab},
		  _gauss_points(gauss_points),
		  _Xi(gauss_points), _Wi(gauss_points)
	{
		if (!std::ranges::is_sorted(xtab)) {
			log(LOG_WARNING, "Grid", "xtab array was not sorted. will sort it and continue");
			std::ranges::sort(_xtab);
		}
		
		switch (grid_fill_type)
		{
			case 1: initGrid(_xtab, nx); break;
			case 2: initGrid2(_xtab, nx); break;
			case 3: initGrid3(_xtab, nx); break;
			default: {
				log(LOG_WARNING, "Grid", "Invalid grid fill type. Found {}, expected 1, 2, or 3.", grid_fill_type);
				log(LOG_WARNING, "Grid", "Will default to 1 (candia-v1 method)");
				initGrid(xtab, nx);
			}
		}
	    
		initGauLeg(0.0, 1.0, _Xi, _Wi);
	}

	Grid::Grid(Grid& other)
		: _points{other._points}, _ntab{other._ntab}, _xtab{other._xtab},
		  _gauss_points{other._gauss_points}, _Xi{other._Xi}, _Wi{other._Wi},
		  _split_interval{other._split_interval}, _Xi2{other._Xi2}, _Wi2{other._Wi2},
		  _gsl_gauleg_table{std::move(other._gsl_gauleg_table)}
	{
		log(LOG_INFO, "Grid", "Copy constructor: previous grid's gsl objects will no longer be valid.");
	}

	void Grid::operator=(Grid& other)
	{
		_points = other._points;
		_ntab = other._ntab;
		_xtab = other._xtab;
		_gauss_points = other._gauss_points;
		_Xi = other._Xi;
		_Wi = other._Wi;
		_split_interval = other._split_interval;
		_Xi2 = other._Xi2;
		_Wi2 = other._Wi2;
		_gsl_gauleg_table = std::move(other._gsl_gauleg_table);
		
		log(LOG_INFO, "Grid", "Copy constructor: previous grid's gsl objects will no longer be valid.");
	}

	void Grid::initGrid(std::vector<double> const& xtab, const uint nx)
	{
		log(LOG_INFO, "Grid", "Using Candia-v1 method.");
		
		const uint xtab_len = xtab.size();
		std::vector<double> Ntab(xtab_len);
	    std::vector<int> ntab(xtab_len);

		double temp = -std::log10(xtab[0]);

		for (uint i=1; i<xtab_len; i++)
			Ntab[i] = (double)(nx-1)*std::log10(xtab[i]/xtab[i-1])/temp;

		ntab[0] = nx-1;

		for (uint i=1; i<xtab_len; i++) {
			ntab[i] =  (int)Ntab[i];
			Ntab[i] -= (double)ntab[i];
			ntab[0] -= ntab[i];
		}

		for (uint i=1; i<xtab_len; i++) {
			if (ntab[i] == 0) {
				ntab[i] =  1;
				Ntab[i] -= 1.0;
				ntab[0] -= 1;
			}
		}

		uint n;
		for ( ; ntab[0]<0; ntab[0]++) {
			n=0;
			for (uint i=1; i<xtab_len; i++) {
				if (ntab[i] != 1) {
					if ((n == 0) || (Ntab[i] <= temp)) {
						n = i;
						temp = Ntab[i];
					}
				}
			}

			ntab[n]--;
			Ntab[n] += 1.0;
		}

		for ( ; ntab[0]>0; ntab[0]--) {
			n=0;
			for (uint i=1; i<xtab_len; i++) {
				if ((n == 0) || (Ntab[i] > temp)) {
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
		for (uint i=0; i<xtab_len-1; i++) {
			lstep = std::log10(xtab[i+1]/xtab[i])/(double)(ntab[i+1]-ntab[i]);

			for (int j=ntab[i]; j<ntab[i+1]; j++)
				_points.at(j) = xtab[i]*std::pow(10.0, lstep*(double)(j-ntab[i]));
		}

		_points.at(nx-1) = 1.0;
		_ntab = ntab;

		return;
	}



	void Grid::initGrid2(std::vector<double> const& xtab, uint nx)
	{
		UNUSED(xtab);
		UNUSED(nx);
		
	    log(LOG_INFO, "Grid", "Using method 2");
		log(LOG_WARNING, "Grid", "Method 2 is unfinished. Prefer method 3 for now.");
		log(LOG_WARNING, "Grid", "This method hard-codes values to compare with distributions.");
		log(LOG_WARNING, "Grid", "Will ignore supplied x-tab and number of grid points.");

	    uint grid_points_per = 500;
		
		double log_min = std::log10(1e-5);
		double log_max = std::log10(0.1);
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = grid_points_per/num_log_intervals;
		
		double lin1_min = 0.1;
		double lin1_max = 0.7;
		double lin2_min = 0.7;
		double lin2_max = 1.0;

		_points.clear();
		for (uint i=0; i<num_log_intervals; ++i) {
			double l0 = log_min + i*dlog;
			double l1 = l0 + dlog;
			for (uint k=0; k<log_interval_size; ++k) {
				double l = l0 + (l1-l0)*k/static_cast<double>(log_interval_size);
				_points.emplace_back(std::pow(10, l));
			}
		}

		for (uint k=0; k<grid_points_per; ++k) {
		    double x = lin1_min + (lin1_max-lin1_min)*k/static_cast<double>(grid_points_per);
			_points.emplace_back(x);
		}
		for (uint k=0; k<grid_points_per; ++k) {
		    double x = lin2_min + (lin2_max-lin2_min)*k/static_cast<double>(grid_points_per-1);
			_points.emplace_back(x);
		}
		
		std::set<double> points_set(_points.begin(), _points.end());
		points_set.insert(xtab.begin(), xtab.end());
		std::vector<double> pivots{0.1, 0.7};
		points_set.insert(pivots.begin(), pivots.end()); // just in case
		_points = std::vector<double>(points_set.begin(), points_set.end());
		

		_ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(_points.begin(), _points.end(), x); it != _points.end()) {
				_ntab.emplace_back(std::distance(_points.begin(), it));
				continue;
			}
		}

		log(LOG_DEBUG, "Grid", "Printing the grid:");
		for (double x : _points)
			log(LOG_DEBUG, "Grid", "  {}", x);
		std::ostringstream ss{};
		std::copy(_ntab.begin(), _ntab.end(), std::ostream_iterator<uint>(ss, ", "));
		log(LOG_DEBUG, "Grid", "Printing ntab array: {}", ss.str());
		std::vector<double> xtabbed_points{};
		std::transform(
			_ntab.begin(), _ntab.end(),
			std::back_insert_iterator(xtabbed_points),
			[&](int x) -> double { return this->_points[x]; });
		ss = {};
		std::copy(xtabbed_points.begin(), xtabbed_points.end(), std::ostream_iterator<double>(ss, ", "));
		log(LOG_DEBUG, "Grid", "Printing xtabbed points to make sure ntab is correct: {}", ss.str());

		initGauLeg(0.0, 1.0, _Xi, _Wi);
	}

	void Grid::initGrid3(std::vector<double> const& xtab, uint nx)
	{
		log(LOG_INFO, "Grid", "Using method 3");
		
		std::vector<double> points{};

		std::vector<double> log_tab{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 0.8, 1.0};
		std::vector<double> log_xtab{log_tab};
		std::ranges::transform(log_xtab, log_xtab.begin(), [](double x) -> double{ return std::log10(x); });
		int num_grid_points_per_bin = nx / xtab.size();

		for (uint i=0; i<log_xtab.size()-1; ++i) {
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
		set.insert(xtab.begin(), xtab.end());
		points = std::vector<double>(set.begin(), set.end());
		
		_points = points;

		// build the ntab array
		_ntab = std::vector<int>{};
		for (const double x : xtab) {
			auto it = std::ranges::lower_bound(_points, x);
			if (it != _points.end() && std::abs(*it - x) < 1e-14)
				_ntab.emplace_back(std::distance(_points.begin(), it));
		}
	}

	void Grid::initGauLeg(double x1, double x2, std::vector<double>& Xi, std::vector<double>& Wi)
	{
		const double eps = 3.0e-11; // relative precision

		// abscissae are symmetric:
		uint n = _gauss_points; // simpler to type
		double N = static_cast<double>(n);
		uint m = (n+1)/2;
		// double x2 = 1.0;
		// double x1 = 0.0;
		double xm = 0.5*(x2+x1);
		double xl = 0.5*(x2-x1);

		for (uint i=1; i<=m; i++) {
			double I = static_cast<double>(i);
			double z = std::cos(PI*(I-0.25)/(N+0.5));

			// default initialize some of these.
			// easy to spot if error occurs
			double z1 = std::numeric_limits<double>::max();
			double pp = std::numeric_limits<double>::max();

			double p1, p2, p3;
			double J;
			do {
				p1 = 1.0;
				p2 = 0.0;

				for (uint j=1; j<=n; j++) {
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
				log(LOG_ERROR, "Grid::initGauLeg()", "Failed to determine gauss-legendre abscissae/weights");

			Xi[i-1] = xm - xl*z;
			Xi[n-i] = xm + xl*z;
			Wi[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp);
			Wi[n-i] = Wi[i-1];
		}
	}

	void Grid::splitConvolution(std::vector<double> const& intervals)
	{
		_split_interval = true;

		auto num_intervals = intervals.size();
		if (num_intervals > 5)
			log(LOG_WARNING, "Grid::splitConvolution()", ">5 convolution intervals is a little excessive.");
		
		_Xi2.reserve(num_intervals);
		_Wi2.reserve(num_intervals);

		if (auto it = std::find_if(intervals.begin(), intervals.end(), [](double x) { return (x < 1e-7) || (x > 1.0); });
			it != std::ranges::end(intervals)) {
			log(LOG_ERROR, "Grid::splitConvolution()", "A provided point ({}) is outside the range [1e-7, 1.0]", *it);
		}

		auto intervals_sorted = intervals;
		std::sort(intervals_sorted.begin(), intervals_sorted.end());
		auto points_per_interval = size()/num_intervals;

		auto it = intervals_sorted.begin() + 1;
		while (it != intervals_sorted.end()) {
			auto prev_x = *(it-1);
			auto new_x = *it;
			gauleg_type Xi(_gauss_points, 0.0), Wi(_gauss_points, 0.0);
			initGauLeg(prev_x, new_x, Xi, Wi);
			_Xi2.emplace_back(Xi);
			_Wi2.emplace_back(Wi);
			++it;
		}
	}

    int Grid::interpFindIdx(double x)
	{
		const static int n = static_cast<int>(size());
		const static int max_k = static_cast<int>(n-2*INTERP_POINTS);

		auto it = std::upper_bound(_points.begin(), _points.end(), x);
		int k = static_cast<int>(it - _points.begin()) - INTERP_POINTS;
		
		if (k < 0)
			k = 0;
		else if (k > max_k)
			k = max_k;

		return k;
	}

	double Grid::interpolate(ArrayGrid& yy, double x)
	{
		const static int n = 2*INTERP_POINTS;
		int ns=0;
		double y, den, dif, dift, ho, hp, w;

		int k = interpFindIdx(x);

		double const* xa = &(_points.data()[k]);
		double const* ya = &(yy.base().data()[k]);
		
		std::array<double, 2*INTERP_POINTS> c{};
		std::array<double, 2*INTERP_POINTS> d{};

		dif = std::abs(x - xa[0]);

		for (int i=0; i<n; i++) {
			if ((dift = std::abs(x - xa[i])) < dif) {
				ns = i;
				dif = dift;
			}
			c[i] = ya[i];
			d[i] = ya[i];
		}

		y = ya[ns--];

		for (int m=1; m<n; m++) {
		    for (int i=0; i<n-m; i++) {
				ho = xa[i] - x;
				hp = xa[i+m] - x;
				w = c[i+1] - d[i];

				den = ho-hp;
				if (std::abs(ho-hp) < 1e-15)
					log(LOG_ERROR, "Grid::interpolate()", "found a denominator equal to 0.0 ({}, {}, {})", x, xa[i], xa[i+m]);

				den = w/den;
				d[i] = hp*den;
				c[i] = ho*den;
			}

			y += (2*ns < (n-1-m) ? c[ns+1] : d[ns--]);
		}

		return y;
	}

	double Grid::convolution(ArrayGrid& A, Expression &E, uint k)
	{
		double x = _points[k];
		double logx =  std::log(x);
		double eplus1 = E.plus(1.0);
		double ed1 = E.delta(1.0);
		double res = (eplus1*std::log1p(-x) + ed1) * A[k];

		auto conv_interval = [&](gauleg_type const& X, gauleg_type const& W) {
			for (uint i=0; i<_gauss_points; i++) {
				double y = X[i];
				double w = W[i];

				double a = std::pow(x, 1.0-y);
				double b = std::pow(x, y);

				double interp1 = interpolate(A, b);
				double interp2 = interpolate(A, a);

				double erega = E.regular(a);
				double eplusb = E.plus(b);

				res -= w*logx*a*erega*interp1;
				res -= w*logx*b*(eplusb*interp2 - eplus1*A[k])/(1.0-b);
			}
		};
		
		if (!_split_interval)
			conv_interval(_Xi, _Wi);
		else {
			for (uint i=0; i<_Xi2.size(); ++i) {
				auto const& X = _Xi2[i];
				auto const& W = _Wi2[i];
				conv_interval(X, W);
			}
		}
	    
		return res;
	}
}
