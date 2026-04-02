#include "Candia-v2/LHAPDFDataFile.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <format>
#include <iterator>

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

		fs::path infofile(dir);
		infofile /= _name;
		infofile.replace_extension(".info");
		log(LOG_DEBUG, "LHAPDFDataFile", "Trying to open file '{}'", infofile.string());
		std::ofstream info_file_stream(infofile);
		if (!info_file_stream)
			log(LOG_ERROR, "LHAPDFDataFile", "Failed to create the datafile with path: '{}'", infofile.string());

		std::string infofile_in_contents = read_file(_infofile_path);
		std::copy(infofile_in_contents.begin(), infofile_in_contents.end(), std::ostreambuf_iterator<char>(info_file_stream));

		fs::path datafile(dir);
		std::string datafile_name = std::format("{}_0000.dat", _name);
		datafile /= datafile_name;
	    log(LOG_DEBUG, "LHAPDFDataFile", "Trying to open file '{}'", datafile.string());
		std::ofstream datafile_stream(datafile);
		if (!datafile_stream)
			log(LOG_ERROR, "LHAPDFDataFile", "Failed to create the datafile with path: '{}'", infofile.string());
		datafile_stream << "test";
	}
}
