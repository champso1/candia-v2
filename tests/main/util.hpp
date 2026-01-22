#ifndef __UTIL_HPP
#define __UTIL_HPP

#include <functional>
#include <string>
#include <cmath>
#include <format>
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

#include "Candia-v2/Common.hpp"
using namespace Candia2;
using dist_type = std::vector<std::vector<double>>;
using xtab_type = std::vector<double>;

static std::string scientificToLatex(double num, int precision, bool benchmark_format)
{
	int exponent = std::floor(std::log10(std::abs(num)));
	double mantissa = num / std::pow(10, exponent);
	if (benchmark_format) {
		return std::vformat("${0: .{1}f}^{{{2:+}}}$",
			std::make_format_args(mantissa, precision, exponent));
	} else {
		int precision_new = precision-1;
		return std::vformat("{0: .{1}f}e${2:+}$",
			std::make_format_args(mantissa, precision_new, exponent));
	}
}
static std::string percentToLatex(double num)
{
	double percent = num*100.0;
	if (percent >= 1000.0) {
		int exponent = std::floor(std::log10(std::abs(num)));
		double mantissa = num / std::pow(10, exponent);
		return std::format("${:.2f}^{{{:+}}}$\\%", mantissa, exponent);
	} else {
    	return std::format("${:.2f}$\\%", num*100.0);
	}
}

static void file_exists(fs::path const& path)
{
	if (!fs::exists(path))
		log(LOG_ERROR, "compare.cpp", "Failed to find file '{}'", path.string());
}

static void print_dist(dist_type const& dist) {
	for (int j=0; j<dist.at(0).size(); ++j) {
		for (int i=0; i<dist.size(); ++i)
			std::cout << dist[i][j] << ' ';
		std::cout << '\n';
	}
	std::cout << std::endl;
}

enum CompareType
{
	ALL_FLAVORS = 0,
	SPECIAL_COMBOS = 1,
	SPECIAL_COMBOS_QM = 2
};
static std::vector<std::string> cols_all_flavors{"g", "xu", "xd", "xs", "xc", "xb", "xub", "xdb", "xsb", "xcb", "xbb"};
static std::vector<std::string> cols_special_combos{"xuv", "xdv", "xL-", "xL+", "xs+", "xc+", "xb+", "xg"};
static std::vector<std::string> cols_special_combos_qm{"xuv", "xdv", "xL-", "xL+", "xs+", "xc+", "xb+", "xg", "xq(-)"};
static std::vector<std::reference_wrapper<const std::vector<std::string>>> cols{
	std::cref(cols_all_flavors),  std::cref(cols_special_combos), std::cref(cols_special_combos_qm)};

static std::pair<xtab_type, dist_type> read_candia_file(fs::path const &path, int size)
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

	while (std::getline(file, line)) {
		iss = std::istringstream(line);
		iss >> temp;
		for (int i=0; i<size; ++i) {
			iss >> temp;
			dists.at(i).push_back(temp);
		}

	}
	dist_type dists_ntabbed(size, std::vector<double>(ntab.size()-1, 0.0));
	for (uint i=0; i<size; ++i) {
		for (uint j=0; j<ntab.size()-1; ++j) {
			uint idx = ntab[j];
			dists_ntabbed[i][j] = dists[i][idx];
		}
	}
	return {xtab, dists_ntabbed};
}

static dist_type fix_dists(dist_type const& dists, int type)
{
	uint ncols = cols[type].get().size();
	switch (type) {
		case 0: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (int k=0; k<dists_fixed.at(0).size(); ++k) {
				for (int j=0; j<dists_fixed.size(); ++j) {
					int idx = j;
					if (j > 5)
						idx += 1;
					dists_fixed.at(j).at(k) = dists[idx][k];
				}
			}
			return dists_fixed;
		} break;
		case 1: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (int k=0; k<dists_fixed.at(0).size(); ++k) {
				dists_fixed.at(0).at(k) = dists[1][k] - dists[6+1][k];
				dists_fixed.at(1).at(k) = dists[2][k] - dists[6+2][k];
				dists_fixed.at(2).at(k) = dists[6+2][k] - dists[6+1][k];
				dists_fixed.at(3).at(k) = 2.0*(dists[2+6][k] + dists[6+1][k]);
				dists_fixed.at(4).at(k) = dists[3][k] + dists[6+3][k];
				dists_fixed.at(5).at(k) = dists[4][k] + dists[6+4][k];
				dists_fixed.at(6).at(k) = dists[5][k] + dists[6+5][k];
				dists_fixed.at(7).at(k) = dists[0][k];
			}
			return dists_fixed;
		} break;
		case 2: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (int k=0; k<dists_fixed.at(0).size(); ++k) {
				dists_fixed.at(0).at(k) = dists[1][k] - dists[6+1][k];
				dists_fixed.at(1).at(k) = dists[2][k] - dists[6+2][k];
				dists_fixed.at(2).at(k) = dists[2+6][k] - dists[6+1][k];
				dists_fixed.at(3).at(k) = 2.0*(dists[2+6][k] + dists[6+1][k]);
				dists_fixed.at(4).at(k) = dists[3][k] + dists[6+3][k];
				dists_fixed.at(5).at(k) = dists[4][k] + dists[6+4][k];
				dists_fixed.at(6).at(k) = dists[5][k] + dists[6+5][k];
				dists_fixed.at(7).at(k) = dists[0][k];
				dists_fixed.at(8).at(k) = 0.0;

				for (uint j=0; j<=6; ++j)
					dists_fixed.at(8).at(k) += dists[j][k] - dists[j+6][k];
			}
			return dists_fixed;
		}
	}
	return {};
}

static std::pair<xtab_type, dist_type> read_other_file(fs::path const &path, int size)
{
	dist_type dists(size, std::vector<double>{});
	std::ifstream file(path);
	file.ignore(1000, '\n');

	xtab_type xtab{};
	std::string line;
	double temp;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		iss >> temp;
		xtab.emplace_back(temp);
		for (int i=0; i<size; ++i) {
			iss >> temp;
			dists.at(i).emplace_back(temp);
		}
	}
	return {xtab, dists};
}

constexpr static char const* TEX_TABLE_DIR{"tex-table"};
constexpr static char const* TEX_TABLE_TEMPLATE{"table-base.txt"};
constexpr static char const* TEX_SUBTABLE_TEMPLATE{"table-sub-base.txt"};
constexpr static char const* TEX_FOOTER_TEMPLATE{"table-footer.txt"};
constexpr static char const* TEX_TABLE_COL_LINE{"    \\multicolumn{1}{c|} {$^COL^$} &\n"};
constexpr static char const* TEX_TABLE_COL_LINE_FINAL{"    \\multicolumn{1}{c||}{$^COL^$} \\\\[0.5mm]"};
constexpr static char const* TEX_TABLE_COL_DEF{"r|"};

static void outputLatexTable(
	std::vector<double> const& xtab,
	dist_type const& diffs, std::string const& filename,
	std::vector<std::string> const& cols, bool use_percentages, bool benchmark_format)
{
	fs::path tex_table_dir = fs::current_path()/TEX_TABLE_DIR;
	fs::path tex_table_base = tex_table_dir/TEX_TABLE_TEMPLATE;
	fs::path tex_subtable = tex_table_dir/TEX_SUBTABLE_TEMPLATE;
    fs::path tex_table_footer = tex_table_dir/TEX_FOOTER_TEMPLATE;
	if (!exists(tex_table_base) || !exists(tex_subtable) || !exists(tex_table_footer))
		log(LOG_ERROR, "util.hpp", "Failed to open the tex template files.");
	
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
	while (pos != std::string::npos) {
		main_table.replace(pos, 6, ncols);
		pos = main_table.find("^COLS^", pos);
	}
	for (std::string const& col : cols | std::ranges::views::take(cols.size()-1)) {
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
	while (pos != std::string::npos) {
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
	if (!fs::exists(latex_build_dir)) {
		if (!fs::create_directory(latex_build_dir))
			log(LOG_ERROR, "util.hpp", "Failed to create latex build directory.");
		log(LOG_INFO, "util.hpp", "'latex' directory created.");
	} else {
	    log(LOG_INFO, "util.hpp", "'latex' directory exists. Continuing.");
	}

	auto format_val = [b=use_percentages,f=benchmark_format](double val) -> std::string {
		if (b)
			return percentToLatex(val);
		else
			return scientificToLatex(val, 4, f);
	};

	log(LOG_INFO, "util.hpp", "Printing table information.");
	for (int i=0; i<diffs.at(0).size(); ++i) {
		double x = xtab.at(i);
		table_text += scientificToLatex(x, 1, benchmark_format) + " & ";
				
		for (uint j=0; j<diffs.size()-1; ++j)
			table_text += format_val(diffs.at(j).at(i)) + " & ";
		table_text += format_val(diffs.back().at(i));
			
		table_text += " \\\\\n";
	}

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
	log(LOG_INFO, "util.hpp", "Cleaning up auxilliary files...");

	fs::path pdf_path(fs::current_path()/fs::path("latex")/fs::path(filename + ".pdf")), new_pdf_path{fs::current_path()};
	fs::copy(pdf_path, new_pdf_path, fs::copy_options::overwrite_existing);
	auto dir_view =
		fs::directory_iterator{fs::current_path()/"latex"}
		| std::ranges::views::filter([&filename](fs::directory_entry const& e) -> bool {
			if (e.path().has_extension() && e.path().extension().string() != ".tex" && e.path().filename().string().starts_with(filename))
				return true;
			return false;
		});
	for (fs::directory_entry const& e : dir_view)
		fs::remove(e.path());
}




#endif // __UTIL_HPP
