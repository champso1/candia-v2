#include <charconv>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <print>
#include <sstream>
#include <string>
#include <ranges>
#include <vector>
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;
using dist_type = std::vector<std::vector<double>>;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

constexpr static char const* TEX_TABLE_DIR{"tex-table"};
constexpr static char const* TEX_TABLE_TEMPLATE{"table-base.txt"};
constexpr static char const* TEX_SUBTABLE_TEMPLATE{"table-sub-base.txt"};
constexpr static char const* TEX_FOOTER_TEMPLATE{"table-footer.txt"};
constexpr static char const* TEX_TABLE_COL_LINE{"    \\multicolumn{1}{c|} {$^COL^$} &\n"};
constexpr static char const* TEX_TABLE_COL_LINE_FINAL{"    \\multicolumn{1}{c||}{$^COL^$} \\\\[0.5mm]"};
constexpr static char const* TEX_TABLE_COL_DEF{"r|"};

static std::vector<double> XTAB{1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};

static std::string scientificToLatex(double num, uint precision)
{
	int exponent = std::floor(std::log10(num));
	double mantissa = num / std::pow(10, exponent);
	return std::format("${:.1f}^{{{}}}", mantissa, exponent);
}
static std::string percentToLatex(double num)
{
    return std::format("{:.2f}\\%", num*100.0);
}

void outputLatexTable(dist_type const& diffs, std::string const& filename, std::vector<std::string> const& cols, std::string const& caption="");


static void usage()
{
	std::println("USAGE: ./compare <candia-file> <other-file> <type>");
	std::println("    <type>: 0=all flavors independently, 1=special combos from benchmark paper");
	exit(EXIT_FAILURE);
}

static void file_exists(fs::path const& path)
{
	if (!fs::exists(path))
	{
		println("[ERROR] compare.cpp: Failed to find file '{}'", path.string());
		exit(EXIT_FAILURE);
	}
}


dist_type read_candia_file(fs::path const &path, int size=13);
dist_type fix_candia_dists(dist_type const& candia, int type);
dist_type read_other_file(fs::path const &path, int size);
dist_type compute_diffs(dist_type const& candia, dist_type const& other);

enum CompareType
{
	ALL_FLAVORS = 0,
	SPECIAL_COMBOS = 1
};

int main(int argc, char *argv[])
{
	if (argc != 4)
		usage();

	fs::path candia_filepath{argv[1]}, other_filepath{argv[2]};
	
	int type;
	std::from_chars(argv[3], argv[3] + 1, type);
	if (type != 0 && type != 1)
	{
		std::println("[ERROR] compare.cpp: Invalid type: {}", type);
		usage();
	}
	std::vector<std::string> cols_all_flavors{"g", "xu", "xd", "xs", "xc", "xb", "xub", "xdb", "xsb", "xcb", "xbb"};
	std::vector<std::string> cols_special_combos{"xuv", "xdv", "xL-", "xL+", "xs+", "xc+", "xb+", "xg"};
	const int ncols = type == 0 ? cols_all_flavors.size() : cols_special_combos.size();
	
    file_exists(candia_filepath);
	file_exists(other_filepath);

	dist_type candia_dists_raw = read_candia_file(candia_filepath);
	dist_type other_dists = read_other_file(other_filepath, ncols);

	if (candia_dists_raw.at(0).size() != other_dists.at(0).size())
	{
		std::println("[ERROR] compare.cpp: Data size mismatch.");
		std::println("[INFO^] compare.cpp: Candia size: {}, other size: {}", candia_dists_raw.at(0).size(), other_dists.at(0).size());
		exit(EXIT_FAILURE);
	}
	
	dist_type candia_dists = fix_candia_dists(candia_dists_raw, type);

	if (candia_dists.size() != other_dists.size())
	{
		std::println("[ERROR] compare.cpp: Dist size mismatch.");
		std::println("[INFO^] compare.cpp: Candia size: {}, other size: {}", candia_dists.size(), other_dists.size());
		exit(EXIT_FAILURE);
	}

	auto diffs = compute_diffs(candia_dists, other_dists);

    
	std::string identifier = candia_filepath.filename().string().substr(0, candia_filepath.filename().string().rfind('.'));
	std::string latex_filename = std::format("diffs-other-{}", identifier);
	outputLatexTable(
		diffs, latex_filename,
		type == 0 ? cols_all_flavors : cols_special_combos);
}



dist_type read_candia_file(fs::path const &path, int size)
{
	dist_type dists(size, std::vector<double>{});
	std::ifstream file(path);
	
	std::vector<double> xtab{};
	std::vector<int> ntab{};
	double temp;
	int temp2;
	std::string line{};

	file.ignore(1000, '\n');
	
	std::getline(file, line);
	std::istringstream iss(line);
	while (iss >> temp)
		xtab.push_back(temp);
	
	std::getline(file, line);
	iss = std::istringstream(line);
	while (iss >> temp2)
		ntab.push_back(temp2);

	while (std::getline(file, line))
	{
		iss = std::istringstream(line);
		iss >> temp;
		for (int i=0; i<size; ++i)
		{
			iss >> temp;
			dists.at(i).push_back(temp);
		}

	}
	dist_type dists_ntabbed(size, std::vector<double>(ntab.size()-1, 0.0));
	for (uint i=0; i<size; ++i)
	{
		for (uint j=0; j<ntab.size()-1; ++j)
		{
			uint idx = ntab[j];
			dists_ntabbed[i][j] = dists[i][idx];
		}
	}
	return dists_ntabbed;
}

dist_type fix_candia_dists(dist_type const& candia, int type)
{
	if (type == 0)
	{
		dist_type candia_dists(11, std::vector<double>(candia.at(0).size(), 0.0));
		for (int k=0; k<candia_dists.at(0).size(); ++k)
		{
			candia_dists.at(0).at(k) =  candia[0][k];
			candia_dists.at(1).at(k) =  candia[1][k];
			candia_dists.at(2).at(k) =  candia[2][k];
			candia_dists.at(3).at(k) =  candia[3][k];
			candia_dists.at(4).at(k) =  candia[4][k];
			candia_dists.at(5).at(k) =  candia[5][k];
			candia_dists.at(6).at(k) =  candia[6+1][k];
			candia_dists.at(7).at(k) =  candia[6+2][k];
			candia_dists.at(8).at(k) =  candia[6+3][k];
			candia_dists.at(9).at(k) =  candia[6+4][k];
			candia_dists.at(10).at(k) = candia[6+5][k];
		}
		return candia_dists;
	}
	else
	{
		dist_type candia_dists(8, std::vector<double>(candia.at(0).size(), 0.0));
		for (int k=0; k<candia_dists.at(0).size(); ++k)
		{
			candia_dists.at(0).at(k) = candia[1][k] - candia[6+1][k];
			candia_dists.at(1).at(k) = candia[2][k] - candia[6+2][k];
			candia_dists.at(2).at(k) = candia[2+6][k] - candia[6+1][k];
			candia_dists.at(3).at(k) = 2.0*(candia[2+6][k] + candia[6+1][k]);
			candia_dists.at(4).at(k) = candia[3][k] + candia[6+3][k];
			candia_dists.at(5).at(k) = candia[4][k] + candia[6+4][k];
			candia_dists.at(6).at(k) = candia[5][k] + candia[6+5][k];
			candia_dists.at(7).at(k) = candia[0][k];
		}
		return candia_dists;
	}
}

dist_type read_other_file(fs::path const &path, int size)
{
	dist_type dists(size, std::vector<double>{});
	std::ifstream file(path);
	file.ignore(1000, '\n');

	std::string line;
	double temp;
	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		iss >> temp;
		for (int i=0; i<size; ++i)
		{
			iss >> temp;
			dists.at(i).emplace_back(temp);
		}
	}
	return dists;
}

dist_type compute_diffs(dist_type const& candia_data, dist_type const& other_data)
{
	auto reldiff =
		[](double candia, double other) -> double {
			return std::abs(candia-other)/other;
		};

	dist_type diffs{candia_data.size(), std::vector<double>(candia_data.at(0).size(), 0.0)};
	for (uint j=0; j<candia_data.size(); ++j)
	{
		for (uint k=0; k<candia_data.at(0).size(); ++k)
		{
			double candia = candia_data.at(j).at(k);
			double other = other_data.at(j).at(k);
			diffs.at(j).at(k) = reldiff(candia, other);
		}
	}
	return diffs;
}

void outputLatexTable(dist_type const& diffs, std::string const& filename, std::vector<std::string> const& cols, std::string const& caption)
{
	fs::path tex_table_dir = fs::current_path()/TEX_TABLE_DIR;
	fs::path tex_table_base = tex_table_dir/TEX_TABLE_TEMPLATE;
	fs::path tex_subtable = tex_table_dir/TEX_SUBTABLE_TEMPLATE;
    fs::path tex_table_footer = tex_table_dir/TEX_FOOTER_TEMPLATE;
	if (!exists(tex_table_base) || !exists(tex_subtable) || !exists(tex_table_footer))
	{
		std::println("[ERROR] compare.cpp: failed to open the tex template files.");
		exit(EXIT_FAILURE);
	}
	
	std::string ncols = std::format("{}", cols.size()+1);
	int pos;
	std::string table_text{};
	
	std::ifstream main_table_s(tex_table_base);
	std::string main_table{std::istreambuf_iterator<char>(main_table_s), std::istreambuf_iterator<char>{}};
	pos = main_table.find("^R^");
	std::string col_def{};
	for (int i=0; i<cols.size(); ++i)
		col_def += TEX_TABLE_COL_DEF;
	main_table.replace(pos, 3, col_def);
	pos = main_table.find("^COLS^");
	while (pos != std::string::npos)
	{
		main_table.replace(pos, 6, ncols);
		pos = main_table.find("^COLS^", pos);
	}
	for (std::string const& col : cols | std::ranges::views::take(cols.size()-1))
	{
		std::string line(TEX_TABLE_COL_LINE);
		pos = line.find("^COL^");
		line.replace(pos, 5, col);
		main_table += line;
	}
	std::string line_final(TEX_TABLE_COL_LINE_FINAL);
	pos = line_final.find("^COL^");
	line_final.replace(pos, 5, cols.back());
	main_table += line_final;
	
	table_text += main_table;
	
    std::ifstream sub_table_s(tex_subtable);
	std::string sub_table{std::istreambuf_iterator<char>(sub_table_s), std::istreambuf_iterator<char>{}};
	pos = sub_table.find("^KR^");
	sub_table.replace(pos, 4, "1.0");
	pos = sub_table.find("^COLS^");
	while (pos != std::string::npos)
	{
		sub_table.replace(pos, 6, ncols);
		pos = sub_table.find("^COLS^", pos);
	}
	std::string amps{};
	for (int i=0; i<cols.size(); ++i)
		amps += " &";
	pos = sub_table.find("^AMPS^");
	sub_table.replace(pos, 6, amps);
	
	table_text += sub_table;
	
	fs::path latex_build_dir = fs::current_path()/"latex";
	if (!fs::exists(latex_build_dir))
	{
		if (!fs::create_directory(latex_build_dir))
		{
			std::println("[ERROR] compare.cpp: Failed to create latex build directory.");
			exit(EXIT_FAILURE);
		}
		std::println("[INFO] compare.cpp: 'latex' directory created.");
	}
	else
	    std::println("[INFO] compare.cpp: 'latex' directory exists. Continuing.");

	std::print("Printing table information...");
	for (int i=0; i<diffs.at(0).size(); ++i)
	{
		double x = XTAB.at(i);
		table_text += scientificToLatex(x, 1) + " & ";
				
		for (uint j=0; j<diffs.size()-1; ++j)
			table_text += percentToLatex(diffs.at(j).at(i)) + " & ";
		table_text += percentToLatex(diffs.back().at(i));
			
		table_text += " \\\\\n";
	}
	std::println("Done.");

	std::ifstream table_footer_s(tex_table_footer);
	std::string table_footer{std::istreambuf_iterator<char>(table_footer_s), std::istreambuf_iterator<char>{}};
	table_text += table_footer;

	std::string title = filename + ".tex";
	fs::path latex_file_path = latex_build_dir/title;
	std::ofstream latex_file(latex_file_path);
	latex_file << table_text;
	latex_file.close();

	std::string command = "pdflatex -interaction=batchmode -output-directory latex " + title;
	system(command.c_str());
	std::println("Cleaning up auxilliary files...");

	fs::path pdf_path(fs::current_path()/fs::path("latex")/fs::path(filename + ".pdf")), new_pdf_path{fs::current_path()};
	fs::copy(pdf_path, new_pdf_path, fs::copy_options::overwrite_existing);
	auto dir_view =
		fs::directory_iterator{fs::current_path()/"latex"}
		| std::ranges::views::filter([&filename](fs::directory_entry const& e) -> bool
		    {
		    	if (e.path().has_extension() && e.path().extension().string() != ".tex" && e.path().filename().string().starts_with(filename))
		    		return true;
		    	return false;
		    });
	for (fs::directory_entry const& e : dir_view)
		fs::remove(e.path());
}
