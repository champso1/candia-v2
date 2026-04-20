#include "Candia-v2/Grid.hpp"

#include <set>
#include <algorithm>
#include <cmath>

namespace Candia2
{
	uint GridFillerLin::fill(std::vector<double>& points)
	{
		// this is to make the intervals less "clean"
		// sometimes, when they are "clean", the linear mapping places points basically right on
		// the xtab points, like 0.3, but off by a delta which is small enough to mess up interpolation
		// if this number is a bit uneven, the hope is that points won't be placed so "cleanly" near
		// xtabbed points, avoiding these errors
		double max = 1.0;

		points.clear();
		for (uint k=0; k<size; ++k) {
		    double x = min + (max-min)*k/static_cast<double>(size);
			points.emplace_back(x);
		}
		
		std::set<double> points_set(points.begin(), points.end());
		points_set.insert(xtab.begin(), xtab.end());
		points = std::vector<double>(points_set.begin(), points_set.end());
		
		ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(points.begin(), points.end(), x); it != points.end()) {
				ntab.emplace_back(std::distance(points.begin(), it));
				continue;
			}
		}

		return points.size();
	}
	std::vector<GridFillerBase::mapping_function_type> GridFillerLin::_mappings
	{
		[](double x, double z) { return std::make_pair(x+(0.9-x)*z, 0.9-x); },
	};
	std::span<GridFillerBase::mapping_function_type> GridFillerLin::_mapping_span(_mappings.begin(), _mappings.end());
	

	uint GridFillerLog::fill(std::vector<double>& points)
	{
		const uint xtab_len = xtab.size();
		std::vector<double> Ntab(xtab_len);
		ntab = std::vector<int>(xtab_len);
		
		points.clear();
		points.resize(size);

		double temp = -std::log10(xtab[0]);

		for (uint i=1; i<xtab_len; i++)
			Ntab[i] = (double)(size-1)*std::log10(xtab[i]/xtab[i-1])/temp;

		ntab[0] = size-1;

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
				points.at(j) = xtab[i]*std::pow(10.0, lstep*(double)(j-ntab[i]));
		}

		points.at(size-1) = 1.0;

		return points.size();
	}
	std::vector<GridFillerBase::mapping_function_type> GridFillerLog::_mappings
	{
		[](double x, double z) {
			auto a = 0.1/x;
			return std::make_pair(x*std::pow(a, z), x*std::pow(a, z)*std::log(a));
		},
	};
	std::span<GridFillerBase::mapping_function_type> GridFillerLog::_mapping_span(_mappings.begin(), _mappings.end());

	
	uint GridFillerLogLin::fill(std::vector<double>& points)
	{
		double log_min = std::log10(min);
		double log_max = std::log10(pivot);
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = log_size/num_log_intervals;

		double lin_min = pivot;
		double lin_max = 1.0;

		points.clear();
		for (uint i=0; i<num_log_intervals; ++i) {
			double l0 = log_min + i*dlog;
			double l1 = l0 + dlog;
			for (uint k=0; k<log_interval_size; ++k) {
				double l = l0 + (l1-l0)*k/static_cast<double>(log_interval_size);
				points.emplace_back(std::pow(10, l));
			}
		}
		
		for (uint k=0; k<lin_size; ++k) {
		    double x = lin_min + (lin_max-lin_min)*k/static_cast<double>(lin_size);
			points.emplace_back(x);
		}
		
		std::set<double> points_set(points.begin(), points.end());
		points_set.insert(xtab.begin(), xtab.end());
		points = std::vector<double>(points_set.begin(), points_set.end());
		
		ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(points.begin(), points.end(), x); it != points.end()) {
				ntab.emplace_back(std::distance(points.begin(), it));
				continue;
			}
		}

		return points.size();
	}

	uint GridFillerLogLinQuad::fill(std::vector<double>& points)
	{
		double log_min = std::log10(min);
		double log_max = std::log10(pivot1);
		uint num_log_intervals = std::round(log_max-log_min);
		double dlog = (log_max-log_min)/static_cast<double>(num_log_intervals);
		uint log_interval_size = log_size/num_log_intervals;

		points.clear();
		for (uint i=0; i<num_log_intervals; ++i) {
			double l0 = log_min + i*dlog;
			double l1 = l0 + dlog;
			for (uint k=0; k<log_interval_size; ++k) {
				double l = l0 + (l1-l0)*k/static_cast<double>(log_interval_size);
				points.emplace_back(std::pow(10, l));
			}
		}

		double lin_min = pivot1;
		double lin_max = pivot2;
		
		for (uint k=0; k<lin_size; ++k) {
		    double x = lin_min + (lin_max-lin_min)*k/static_cast<double>(lin_size);
			points.emplace_back(x);
		}

		double quad_min = pivot2;
		double quad_max = 1.0;

		for (uint k=0; k<quad_size; ++k) {
			double f = 1.0 - k/static_cast<double>(quad_size);
		    double x = quad_max - (quad_max-quad_min)*f*f;
			points.emplace_back(x);
		}
		
		std::set<double> points_set(points.begin(), points.end());
		points_set.insert(xtab.begin(), xtab.end());
		points = std::vector<double>(points_set.begin(), points_set.end());
		
		ntab.clear();
		for (double x : xtab) {
			if (auto it = std::find(points.begin(), points.end(), x); it != points.end()) {
				ntab.emplace_back(std::distance(points.begin(), it));
				continue;
			}
		}
		
		return points.size();
	}
}
