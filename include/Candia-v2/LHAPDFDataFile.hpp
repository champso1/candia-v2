/**
 *  @file LHAPDFDataFile.hpp
 *  @brief Contains the @a LHAPDFDataFile class which handles creating output datafiles consistent with LHAPDF
 */

#ifndef __LHAPDF_DATAFILE_HPP
#define __LHAPDF_DATAFILE_HPP

#include "Candia-v2/Common.hpp"

#include <filesystem>

namespace Candia2
{
	struct InfoFileInParams final
	{
		std::string desc{"Evolution from Candia-v2"};
		std::string authors{"C. Hampson"};
		std::string reference{"arXiv:2512.22667"};
		std::string format{"lhagrid1"};
		int data_version{1};
		int num_members{1};
		int setidx{-1};
		std::vector<int> flavors{-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};
		int order_qcd{2};
		std::string flavor_scheme{"variable"};
		int num_flavors{5};
		std::string error_type{"hessian"};
		int error_conf_level{90};
		double xmin{1e-5};
		double xmax{1.0};
		double qmin{std::numbers::sqrt2};
		double qmax{100.0};
		double mz{MZ};
		double mup{0.0};
		double mdown{0.0};
		double mstrange{std::numbers::sqrt2};
		double mcharm{std::numbers::sqrt2};
	    double mbottom{4.5};
		double mtop{175.0};
		double as_mz{0.118};
		int as_order_qcd{2};
		std::string as_type{"ipol"};
		std::vector<double> as_qs{};
		std::vector<double> as_vals{};
	};
	
	/**
	 *  @brief class to handle creating output datafiles consistent with LHAPDF
	 */
	class LHAPDFDataFile final
	{
		std::string _name;
		std::filesystem::path _infofile_in_path;
	public:
		LHAPDFDataFile(std::string const& name, std::filesystem::path const& infofile_in_path)
			: _name{name}, _infofile_in_path{infofile_in_path} {}

		void setParams(InfoFileInParams const& params);
		void write();
	};
}

#endif // __LHAPDF_DATAFILE_HPP
