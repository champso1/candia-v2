#pragma once

#include <functional>
#include <string>
#include <filesystem>
namespace fs = std::filesystem;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

using dist_type = std::vector<std::vector<double>>;
using xtab_type = std::vector<double>;

std::string percentToLatex(double percent);
std::string scientificToLatex(double num, int precision, bool benchmark_format);
std::string percentToLatex2(double percent);

void file_exists(fs::path const& path);
void print_dist(dist_type const& dist);

enum CompareType
{
	ALL_FLAVORS = 0,
	SPECIAL_COMBOS = 1,
	SPECIAL_COMBOS_QM = 2,
	SPECIAL_COMBOES_NS_AND_S = 3,
	SPECIAL_COMBOES_FFNS = 4,
};

inline std::vector<std::string> cols_all_flavors{"g", "xu", "xd", "xs", "xc", "xb", "x\\bar{u}", "x\\bar{d}", "x\\bar{s}", "x\\bar{c}", "x\\bar{b}"};
inline std::vector<std::string> cols_special_combos{"xuv", "xdv", "xL_-", "xL_+", "xs_+", "xc_+", "xb_+", "xg"};
inline std::vector<std::string> cols_special_combos_qm{"xu_v", "xd_v", "x\\bar{u}", "x\\bar{d}", "xL_+", "xs_+", "xc_+", "xb_+", "xg"};
inline std::vector<std::string> cols_special_combos_ns_and_s{"xu_+", "xc_+", "xb_+", "xq_{NS,1d}^{(+)}", "xq_{NS,1c}^{(+)}", "xq_{NS,1b}^{(+)}", "xq^{(-)}", "xq^{(+)}", "xg"};
inline std::vector<std::string> cols_special_combos_ffns{"xuv", "xdv", "xL-", "xL+", "xs-", "xs+", "xc+", "xg"};
inline std::vector<std::reference_wrapper<const std::vector<std::string>>> cols{
	std::cref(cols_all_flavors),
	std::cref(cols_special_combos),
	std::cref(cols_special_combos_qm),
	std::cref(cols_special_combos_ns_and_s),
	std::cref(cols_special_combos_ffns),
	std::cref(cols_special_combos_ns_and_s),
};

void print_compare_types(std::string_view filename);

std::pair<xtab_type, dist_type> read_candia_file(fs::path const &path, uint size);
dist_type fix_dists(dist_type const& dists, int type);
std::pair<xtab_type, dist_type> read_other_file(fs::path const &path, uint size);

constexpr inline char const* TEX_TABLE_DIR{"tex-table"};
constexpr inline char const* TEX_TABLE_TEMPLATE{"table-base.txt"};
constexpr inline char const* TEX_SUBTABLE_TEMPLATE{"table-sub-base.txt"};
constexpr inline char const* TEX_FOOTER_TEMPLATE{"table-footer.txt"};
constexpr inline char const* TEX_TABLE_COL_LINE{"    \\multicolumn{1}{c|} {$^COL^$} &\n"};
constexpr inline char const* TEX_TABLE_COL_LINE_FINAL{"    \\multicolumn{1}{c||}{$^COL^$} \\\\[0.5mm]"};
constexpr inline char const* TEX_TABLE_COL_DEF{"r|"};

void outputLatexTable(
    xtab_type const& xtab,
	dist_type const& diffs, std::string const& filename,
	std::vector<std::string> const& cols, int format, bool benchmark_format);

