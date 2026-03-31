#include "Candia-v2/LHAPDFDataFile.hpp"

#include <filesystem>
#include <fstream>
#include <format>

namespace Candia2
{
	namespace fs = std::filesystem;
	void LHAPDFDataFile::write()
	{
		fs::path dir(_name);
		if (fs::exists(dir) && fs::is_directory(dir))
			log(LOG_ERROR, "LHAPDFDataFile", "Directory '{}' already exists", dir.string());

		if (!fs::create_directory(dir))
			log(LOG_ERROR, "LHAPDFDataFile", "Failed to create directory '{}'.", dir.string());

		fs::path info_file(dir);
		info_file /= _name;
		info_file.replace_extension(".info");
		log(LOG_DEBUG, "LHAPDFDataFile", "Trying to open file '{}'", info_file.string());
		std::ofstream info_file_stream(info_file);
		if (!info_file_stream)
			log(LOG_ERROR, "LHAPDFDataFile", "Failed to create the datafile with path: '{}'", info_file.string());
		info_file_stream << "test";

		fs::path data_file(dir);
		std::string data_file_name = std::format("{}_0000.dat", _name);
		data_file /= data_file_name;
	    log(LOG_DEBUG, "LHAPDFDataFile", "Trying to open file '{}'", data_file.string());
		std::ofstream data_file_stream(data_file);
		if (!data_file_stream)
			log(LOG_ERROR, "LHAPDFDataFile", "Failed to create the datafile with path: '{}'", info_file.string());
		data_file_stream << "test";
		
	}
}
