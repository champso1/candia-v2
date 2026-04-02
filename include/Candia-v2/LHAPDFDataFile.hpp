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
	/**
	 *  @brief set of some params put into the info file
	 */
	struct InfoFileInParams final
	{
		int order{};
		std::vector<double> as_qs{};
		std::vector<double> as_vals{};
	};
	
	/**
	 *  @brief class to handle creating output datafiles consistent with LHAPDF
	 */
	class LHAPDFDataFile final
	{
		std::string _name;
		std::filesystem::path _infofile_path;
		InfoFileInParams _params;
		
	public:
		LHAPDFDataFile(
			std::string const& name,
			std::filesystem::path const& infofile_path,
			InfoFileInParams&& params)
			: _name{name}, _infofile_path{infofile_path}, _params{params} {}

		void write();
	};
}

#endif // __LHAPDF_DATAFILE_HPP
