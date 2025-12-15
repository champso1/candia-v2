#ifndef __GRID_HPP
#define __GRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <vector>

namespace Candia2
{
	class ArrayGrid;
	class Grid final
	{
	public:
		using grid_type = std::vector<double>;
		using gauleg_type = std::vector<double>;
		using ntab_type = std::vector<int>;

	private:
		grid_type _points{};
		ntab_type _ntab;     //!< stores indices for the tabulated grid points
		
		gauleg_type _Xi{};
		gauleg_type _Wi{};
	public:
		Grid() = delete;
		Grid(grid_type const& xtab, uint nx, int grid_fill_type=1);
		~Grid() = default;


		inline grid_type const& points() const { return _points; }
		inline double at(uint idx) const { return _points.at(idx); };
		inline double operator[](uint idx) const { return _points.at(idx); }

		inline gauleg_type const& abscissae() const { return _Xi; }
		inline double abscissae(uint idx) const { return _Xi[idx]; }
		inline gauleg_type const& weights() const { return _Wi; }
		inline double weights(uint idx) const { return _Wi[idx]; }

		inline uint size() const { return _points.size(); }

		inline ntab_type const& ntab() const { return _ntab; }
		inline int const& ntab(uint idx) const { return _ntab.at(idx); }
		inline ntab_type& ntab() { return _ntab; }
		inline int& ntab(uint idx) { return _ntab.at(idx); }

		inline grid_type::const_iterator begin() const { return _points.begin(); }
		inline grid_type::iterator begin() { return _points.begin(); }
		inline grid_type::const_iterator end() const { return _points.end(); }
		inline grid_type::iterator end() { return _points.end(); }


		uint interpFindIdx(double x) const;
		double interpolate(grid_type const& y, double x);
		double interpolate(ArrayGrid& y, double x);
		double convolution(grid_type const& A, Expression &E, uint k);
		double convolution(ArrayGrid& A, Expression &E, uint k);
		
	private:
		void initGrid(grid_type const& xtab, uint nx);
		void initGrid2(grid_type const& xtab, uint nx);
		void initGrid3(grid_type const& xtab, uint nx);
		void initGauLeg(double x1, double x2, std::vector<double> & Xi, std::vector<double> & Wi);
	};
}

#endif // __GRID_HPP
