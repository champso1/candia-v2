/** @file
 *
 *  Solver for harmonic polylogarithms.
 */

#ifndef __HPLOG_HPP
#define __HPLOG_HPP

#include "Candia-v2/Common.hpp"

#include <iostream>
#include <fstream>
#include <string>


namespace Candia2
{

	class HPLog
	{
	private:
		uint _W; //!< weight with which to calculate

		std::string _file_prefix; ///!< directory to find weight/coeff files

		/** @name arrays for coeffs/weights
		 */
		///@{
		std::vector<double> _DA2, _DA3, _DA4, _DA5, _DA6, _DA7, _DA8;
		const uint LEN[7] =
			{
				360, 1600, 5040, 17280, 51040, 162240, 486000
			};

		std::vector<int> _ii15, _ii25, _ii35, _ii45, _ii55;
		std::vector<int> _ii16, _ii26, _ii36, _ii46, _ii56, _ii66;
		std::vector<int> _ii17, _ii27, _ii37, _ii47, _ii57, _ii67, _ii77;
		std::vector<int> _ii18, _ii28, _ii38, _ii48, _ii58, _ii68, _ii78, _ii88;
		///@}

		/** @brief loads coefficients from file
		 */
		void LoadCoeffs();

		/** @brief loads weights from file
		 */
		void LoadWeights();

		/** @name coefficient checking functions
		 */
		///@{
		int Check_B2(int i1, int i2);
		int Check_B3(int i1, int i2, int i3);
		int Check_B4(int i1, int i2, int i3, int i4);
		int Check_B5(int i1, int i2, int i3, int i4,
					 int i5);
		int Check_B6(int i1, int i2, int i3, int i4,
					 int i5, int i6);
		int Check_B7(int i1, int i2, int i3, int i4,
					 int i5, int i6, int i7);
		int Check_B8(int i1, int i2, int i3, int i4,
					 int i5, int i6, int i7, int i8);
		///@}

		/** @brief helper function to read coeff/weight files
		 */
		template <typename T>
		requires Arithmetic<T>
		std::vector<T> ReadFile(std::string const& filepath, const uint len);

	public:
		/** @brief small enum to indicate what mappings need to take place
		*/
		enum MappingType
		{
			None = 0, //!< no mappings required
		    Negative = 1 << 0, //!< if x is negative, maps it to positive
			GRT1 = 1 << 1, //!< if x is greater than one, maps it to [0,1]
			GRTSqrt2Minus1 = 1 << 2, //!< if greater than sqrt(2)-1, maps it to [0,sqrt(2)-1]
		};

	private:

		/** @brief determines what mappings to do for the given x value
		 */
		uint Check_X(double x);

		/** @name handles mappings if x is out of range
		 */
		///@{
		double H1Mappings(uint mappings, int i1, double x);
		double H2Mappings(uint mappings, int i1, int i2, double x);
		double H3Mappings(uint mappings, int i1, int i2, int i3, double x);
		///@}

		
	public:
		/** @name Constructors/destructors
		 */
		///@{
		HPLog(const uint W, std::string const& file_prefix=".");
		~HPLog() = default;
		///@}

		/** @name actual hplog calclation functions
		 */
		///@{
		double H1(int i1, double x);
		double H2(int i1, int i2, double x);
		double H3(int i1, int i2, int i3, double x);
		double H4(int i1, int i2, int i3, int i4, double x);
		double H5(int i1, int i2, int i3, int i4,
				  int i5, double x);
		double H6(int i1, int i2, int I3, int i4,
				  int i5, int i6, double x);
		double H7(int i1, int i2, int i3, int i4,
				  int i5, int i6, int i7, double x);
		double H8(int i1, int i2, int i3, int i4,
				  int i5, int i6, int i7, int i8,
				  double x);
		///@}
	};


	


	

	// as a templated function, we must place the definition
	// here in the header file
	template <typename T>
	requires Arithmetic<T>
	std::vector<T> HPLog::ReadFile(std::string const& filepath, const uint len)
	{
		std::ifstream f(_file_prefix + '/' + filepath);
		if (!f)
		{
			std::cerr << "[HPLOG] Could not open '" << filepath << "'\n";
			exit(1);
		}

		std::vector<T> vec;
	
		std::string temp;
		uint count = 0;
		while (!f.eof() && count < len) {
			count++;
			f >> temp;
			vec.emplace_back(static_cast<T>(std::stold(temp)));
		}
		f.close();

		return vec;
	}
	
};


#endif // __HPLOG_HPP
