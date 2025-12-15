#ifndef __HPLOG_HPP
#define __HPLOG_HPP

#include "Candia-v2/Common.hpp"

#include <cstdlib>
#include <print>
#include <fstream>
#include <string>
#include <filesystem>

namespace Candia2
{
	class HPLog
	{
	private:
		uint _W; //!< weight with which to calculate
		std::string _file_prefix; //!< directory to find weight/coeff files

		std::vector<double> _DA2, _DA3, _DA4, _DA5, _DA6, _DA7, _DA8;
		const uint LEN[7]{
			360, 1600, 5040, 17280, 51040, 162240, 486000};

		std::vector<int> _ii15, _ii25, _ii35, _ii45, _ii55;
		std::vector<int> _ii16, _ii26, _ii36, _ii46, _ii56, _ii66;
		std::vector<int> _ii17, _ii27, _ii37, _ii47, _ii57, _ii67, _ii77;
		std::vector<int> _ii18, _ii28, _ii38, _ii48, _ii58, _ii68, _ii78, _ii88;


		void loadCoeffs();
		void loadWeights();

		int checkB2(int i1, int i2);
		int checkB3(int i1, int i2, int i3);
		int checkB4(int i1, int i2, int i3, int i4);
		int checkB5(int i1, int i2, int i3, int i4,
					 int i5);
		int checkB6(int i1, int i2, int i3, int i4,
					 int i5, int i6);
		int checkB7(int i1, int i2, int i3, int i4,
					 int i5, int i6, int i7);
		int checkB8(int i1, int i2, int i3, int i4,
					 int i5, int i6, int i7, int i8);

		template <Arithmetic T>
		std::vector<T> readFile(std::string const& filename, uint len);
	public:
		enum MappingType
		{
			None = 0,                //!< no mappings required
		    Negative = 1 << 0,       //!< if x is negative, maps it to positive
			GRT1 = 1 << 1,           //!< if x is greater than one, maps it to [0,1]
			GRTSqrt2Minus1 = 1 << 2, //!< if greater than sqrt(2)-1, maps it to [0,sqrt(2)-1]
		};

	private:

		uint checkX(double x);

		double H1Mappings(uint mappings, int i1, double x);
		double H2Mappings(uint mappings, int i1, int i2, double x);
		double H3Mappings(uint mappings, int i1, int i2, int i3, double x);
		
	public:
		HPLog(const uint W, std::string const& file_prefix=".");
		~HPLog() = default;

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
	};


	template <Arithmetic T>
	std::vector<T> HPLog::readFile(std::string const& filename, uint len)
	{
		namespace fs = std::filesystem;

		fs::path filepath = fs::path{_file_prefix} / fs::path{filename};
		if (!fs::exists(filepath))
		{
			std::println("[HPLOG: ERROR] File '{}' could not be found.", filepath.string());
			exit(EXIT_FAILURE);
		}

		std::ifstream f(filepath);
		std::vector<T> vec;
		std::string temp;
		uint count = 0;
		while (!f.eof() && count < len) {
			count++;
			f >> temp;
			vec.emplace_back(static_cast<T>(std::stold(temp)));
		}
		
		return vec;
	}
	
}; // namespace Candia2


#endif // __HPLOG_HPP
