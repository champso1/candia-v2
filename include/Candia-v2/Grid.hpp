#ifndef __GRID_HPP
#define __GRID_HPP

#include "Candia-v2/Common.hpp"
#include "Candia-v2/Expression.hpp"

#include <memory>
#include <vector>

namespace Candia2
{
	class Grid final
	{
		using grid_type = std::vector<double>;
		using gauleg_type = std::vector<double>;
		using ntab_type = std::vector<int>;
		
		grid_type _points{};
		ntab_type _ntab;     //!< stores indices for the tabulated grid points
		
		gauleg_type _Xi{};
		gauleg_type _Wi{};

		gauleg_type _Xi1{}, _Wi1{};
		gauleg_type _Xi2{}, _Wi2{};
		gauleg_type _Xi3{}, _Wi3{};
	public:

		Grid() = delete;
		Grid(std::vector<double> const& xtab, uint nx, int grid_fill_type=1);
		~Grid() = default;


		inline grid_type const& Points() const { return _points; }
		inline double At(const uint idx) const { return _points.at(idx); };
		inline double operator[](const uint idx) const { return _points.at(idx); }

		inline gauleg_type const& Abscissae() const { return _Xi; }
		inline double Abscissae(uint idx) const { return _Xi[idx]; }
		inline gauleg_type const& Weights() const { return _Wi; }
		inline double Weights(uint idx) const { return _Wi[idx]; }

		inline uint Size() const { return _points.size(); }

		inline ntab_type const& Ntab() const { return _ntab; }
		inline int const& Ntab(uint idx) const { return _ntab.at(idx); }
		inline ntab_type& Ntab() { return _ntab; }
		inline int& Ntab(uint idx) { return _ntab.at(idx); }

		inline grid_type::const_iterator begin() const { return _points.begin(); }
		inline grid_type::iterator begin() { return _points.begin(); }
		inline grid_type::const_iterator end() const { return _points.end(); }
		inline grid_type::iterator end() { return _points.end(); }

		inline gauleg_type const& Abscissae(int idx) const
		{
			switch (idx)
			{
				case 1: return _Xi1;
				case 2: return _Xi2;
				case 3: return _Xi3;
			}
			return _Xi;
		}
		inline gauleg_type const& Weights(int idx) const
		{
			switch (idx)
			{
				case 1: return _Wi1;
				case 2: return _Wi2;
				case 3: return _Wi3;
			}
			return _Wi;
		}


		double Interpolate(grid_type const& y, const double x) const;
		double Convolution(
			grid_type const& A,
			std::shared_ptr<Expression> P,
			uint k);

	private:
		void InitGrid(grid_type const& xtab, const uint nx);
		void InitGrid2(grid_type const& xtab, const uint nx);
		void InitGrid3(grid_type const& xtab, const uint nx);
		void InitGauLeg(double x1, double x2, std::vector<double> & Xi, std::vector<double> & Wi);


	public:
		uint InterpFindIdx(double x) const;
	};
}

#endif // __GRID_HPP
