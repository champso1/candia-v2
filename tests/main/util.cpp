#include "util.hpp"

#include <cmath>

std::string percentToLatex(double percent)
{
	return std::format("${:.2f}$\\%", percent);
}

std::string scientificToLatex(double num, int precision, bool benchmark_format)
{
	int exponent = std::floor(std::log10(std::abs(num)));
	if (std::abs(exponent) > 20)
		return std::format("???");
	double mantissa = num / std::pow(10.0, static_cast<double>(exponent));
	if (benchmark_format) {
		return std::vformat("${0: .{1}f}^{{{2:+}}}$",
			std::make_format_args(mantissa, precision, exponent));
	} else {
		int precision_new = precision-1;
		return std::vformat("{0: .{1}f}e${2:+}$",
			std::make_format_args(mantissa, precision_new, exponent));
	}
}

std::string percentToLatex2(double percent)
{
	if (percent >= 1e-8) {
		int exponent = std::floor(std::log10(std::abs(percent)));
		double mantissa = percent / std::pow(10.0, static_cast<double>(exponent));
		return std::format("{:.2f}e{:+}\\%", mantissa, exponent);
	} else {
		return std::format("0.00\\%");
	}
}

void file_exists(fs::path const& path)
{
	if (!fs::exists(path))
		log(LOG_ERROR, "compare.cpp", "Failed to find file '{}'", path.string());
}

void print_dist(dist_type const& dist) {
	for (uint j=0; j<dist.at(0).size(); ++j) {
		for (uint i=0; i<dist.size(); ++i)
			std::cout << dist[i][j] << ' ';
		std::cout << '\n';
	}
	std::cout << std::endl;
}


void print_compare_types(std::string_view filename)
{
	log(LOG_INFO, filename, "    <type>:      0: all flavors independently");
	log(LOG_INFO, filename, "                 1: special combos from benchmark paper");
	log(LOG_INFO, filename, "                 2: special combos from benchmark paper to avoid cancellations");
	log(LOG_INFO, filename, "                 3: specific non-singlet and singlet distributions");
	log(LOG_INFO, filename, "                 4: ffns distributions");
	log(LOG_INFO, filename, "                 5: same as (3), but using the distributions directly");
	log(LOG_INFO, filename, "                    rather than reconstructing from the individual quark PDFs");
	log(LOG_INFO, filename, "                    requires that the other datafile contains exactly the same info as candia");
}

ReadCandiaFileResult read_candia_file(fs::path const &path, uint size)
{
	dist_type dists(size, std::vector<double>{});
	std::ifstream file(path);
	
	std::vector<double> xtab{};
	std::vector<int> ntab{};
	std::vector<double> grid_points{};
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

	std::getline(file, line);
	iss = std::istringstream(line);
	while (iss >> temp)
		grid_points.push_back(temp);

	while (std::getline(file, line)) {
		iss = std::istringstream(line);
		iss >> temp;
		for (uint i=0; i<size; ++i) {
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
	return {xtab_type(xtab.begin(), xtab.end()-1), dists, dists_ntabbed, grid_points};
}

dist_type fix_dists(dist_type const& dists, int type)
{
	uint ncols = cols[type].get().size();
	switch (type) {
		case 0: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));

			std::vector<int> js{0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11};
			auto tie =
				std::views::iota(0)
				| std::views::take(js.size())
				| std::views::transform([&](int i){ return std::make_pair(i, js[i]); });
		    
			for (uint k=0; k<dists_fixed.at(0).size(); ++k) {
				for (auto [i, j] : tie) {
					dists_fixed.at(i).at(k) = dists[j][k];
				}
			}
			return dists_fixed;
		} break;
		case 1: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (uint k=0; k<dists_fixed.at(0).size(); ++k) {
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
			for (uint k=0; k<dists_fixed.at(0).size(); ++k) {
				dists_fixed.at(0).at(k) = dists[1][k] - dists[6+1][k];
				dists_fixed.at(1).at(k) = dists[2][k] - dists[6+2][k];
				dists_fixed.at(2).at(k) = dists[1+6][k];
				dists_fixed.at(3).at(k) = dists[2+6][k];
				dists_fixed.at(4).at(k) = 2.0*(dists[2+6][k] + dists[6+1][k]);
				dists_fixed.at(5).at(k) = dists[3][k] + dists[3+6][k];
				dists_fixed.at(6).at(k) = dists[4][k] + dists[4+6][k];
				dists_fixed.at(7).at(k) = dists[5][k] + dists[5+6][k];
				dists_fixed.at(8).at(k) = dists[0][k];
			}
			return dists_fixed;
		}; break;
		case 3: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (uint k=0; k<dists_fixed.at(0).size(); ++k) {
				double qpu = dists[1][k] + dists[6+1][k];
				double qmu = dists[1][k] - dists[6+1][k];
				
				dists_fixed.at(0).at(k) = qmu - (dists[2][k] - dists[6+2][k]);
				dists_fixed.at(1).at(k) = qmu - (dists[4][k] - dists[6+4][k]);
				dists_fixed.at(2).at(k) = qmu - (dists[5][k] - dists[6+5][k]);
				dists_fixed.at(3).at(k) = qpu - (dists[2][k] + dists[6+2][k]);
				dists_fixed.at(4).at(k) = qpu - (dists[4][k] + dists[6+4][k]);
				dists_fixed.at(5).at(k) = qpu - (dists[5][k] + dists[6+5][k]);

				dists_fixed.at(6).at(k) = 0.0;
				dists_fixed.at(7).at(k) = 0.0;
				dists_fixed.at(8).at(k) = dists[0][k];

				for (uint j=1; j<=6; ++j) {
					dists_fixed.at(6).at(k) += dists[j][k] - dists[j+6][k];
					dists_fixed.at(7).at(k) += dists[j][k] + dists[j+6][k];
				}
			}
			return dists_fixed;
		}; break;
		case 4: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (uint k=0; k<dists_fixed.at(0).size(); ++k) {
				dists_fixed.at(0).at(k) = dists[1][k] - dists[6+1][k];
				dists_fixed.at(1).at(k) = dists[2][k] - dists[6+2][k];
				dists_fixed.at(2).at(k) = dists[6+2][k] - dists[6+1][k];
				dists_fixed.at(3).at(k) = 2.0*(dists[2+6][k] + dists[6+1][k]);
				dists_fixed.at(4).at(k) = dists[3][k] - dists[6+3][k];
				dists_fixed.at(5).at(k) = dists[3][k] + dists[6+3][k];
				dists_fixed.at(6).at(k) = dists[4][k] + dists[6+3][k];
				dists_fixed.at(7).at(k) = dists[0][k];
			}
			return dists_fixed;
		} break;
		case 5: {
			dist_type dists_fixed(ncols, std::vector<double>(dists.at(0).size(), 0.0));
			for (uint k=0; k<dists_fixed.at(0).size(); ++k) {
			    dists_fixed.at(0).at(k) = dists[19][k];
				dists_fixed.at(1).at(k) = dists[22][k];
				dists_fixed.at(2).at(k) = dists[23][k];
				
				dists_fixed.at(3).at(k) = dists[32][k];
				dists_fixed.at(4).at(k) = dists[34][k];
				dists_fixed.at(5).at(k) = dists[35][k];

				dists_fixed.at(6).at(k) = dists[25][k];
				dists_fixed.at(7).at(k) = dists[31][k];
				dists_fixed.at(8).at(k) = dists[0][k];
			}
			return dists_fixed;
		}; break;
	}
	return dists;
}

std::pair<xtab_type, dist_type> read_other_file(fs::path const &path, uint size)
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
		for (uint i=0; i<size; ++i) {
			iss >> temp;
			dists.at(i).emplace_back(temp);
		}
	}
	return {xtab, dists};
}

void outputLatexTable(
    xtab_type const& xtab,
	dist_type const& diffs, std::string const& filename,
	std::vector<std::string> const& cols, int format, bool benchmark_format)
{
	fs::path tex_table_dir = fs::current_path()/TEX_TABLE_DIR;
	fs::path tex_table_base = tex_table_dir/TEX_TABLE_TEMPLATE;
	fs::path tex_subtable = tex_table_dir/TEX_SUBTABLE_TEMPLATE;
    fs::path tex_table_footer = tex_table_dir/TEX_FOOTER_TEMPLATE;
	if (!exists(tex_table_base) || !exists(tex_subtable) || !exists(tex_table_footer))
		log(LOG_ERROR, "util.hpp", "Failed to open the tex template files.");
	
	std::string ncols = std::format("{}", cols.size()+1);
	std::string::size_type pos;
	std::string table_text{};
	
	std::ifstream main_table_s(tex_table_base);
	std::string main_table{std::istreambuf_iterator<char>(main_table_s), std::istreambuf_iterator<char>{}};
	pos = main_table.find("^R^");
	std::string col_def{};
	for (uint i=0; i<cols.size(); ++i)
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
	for (uint i=0; i<cols.size(); ++i)
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

	auto format_val = [t=format,f=benchmark_format](double val) -> std::string {
		switch (t) {
			case 0: return percentToLatex(val); break;
			case 1: return scientificToLatex(val, 4, f); break;
			case 2: return percentToLatex2(val); break;
		}
		throw std::runtime_error("unreachable");
	};

	log(LOG_INFO, "util.hpp", "Printing table information.");
	for (uint i=0; i<diffs.at(0).size(); ++i) {
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

	auto output_dir = (format == 0 || format == 2) ? "diffs" : "tables";
	fs::path output_dir_path = fs::current_path()/fs::path(output_dir);
	if (!fs::exists(output_dir_path)) {
		if (!fs::create_directory(output_dir_path))
			log(LOG_ERROR, "util.hpp", "Failed to create output directory '{}'.", output_dir_path.string());
	}
	fs::path pdf_path(fs::current_path()/fs::path("latex")/fs::path(filename + ".pdf"));
	log(LOG_INFO, "util.hpp", "{}  ->  {}", pdf_path.string(), output_dir_path.string());
	fs::copy(pdf_path, output_dir_path, fs::copy_options::overwrite_existing);
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

void outputLatexTableScaleRatio(
    xtab_type const& xtab,
	std::array<std::reference_wrapper<dist_type>, 3> const& arrs,
	std::string const& filename,
	std::vector<std::string> const& cols, int format, bool benchmark_format)
{
	fs::path tex_table_dir = fs::current_path()/TEX_TABLE_DIR;
	
	fs::path tex_table_base = tex_table_dir/TEX_TABLE_TEMPLATE;
	fs::path tex_subtable = tex_table_dir/TEX_SUBTABLE_TEMPLATE;
    fs::path tex_table_footer = tex_table_dir/TEX_FOOTER_TEMPLATE;
	if (!exists(tex_table_base) || !exists(tex_subtable) || !exists(tex_table_footer))
		log(LOG_ERROR, "util.hpp", "Failed to open the tex template files.");
	
	std::string ncols = std::format("{}", cols.size()+1);
	std::string::size_type pos;
	std::string table_text{};
	
	std::string main_table = read_file(tex_table_base);
	pos = main_table.find("^R^");
	std::string col_def{};
	for (uint i=0; i<cols.size(); ++i)
		col_def += TEX_TABLE_COL_DEF;
	main_table.replace(pos, 3, col_def);
	pos = main_table.find("^COLS^");
	while (pos != std::string::npos) {
		main_table.replace(pos, 6, ncols);
		pos = main_table.find("^COLS^", pos);
	}
	for (std::string const& col : cols | std::views::take(cols.size()-1)) {
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

	constexpr auto get_mu_str = [&](double mu) -> std::string {
		if (mu == 0.5)
			return std::string("0.5");
		else if (mu == 2.0)
			return std::string("2.0");
		else
			return std::string("1.0");
	};

	auto format_val = [t=format,f=benchmark_format](double val) -> std::string {
		switch (t) {
			case 0: return percentToLatex(val); break;
			case 1: return scientificToLatex(val, 4, f); break;
			case 2: return percentToLatex2(val); break;
		}
		throw std::runtime_error("unreachable");
	};
	
	auto do_subtable = [&](double mu, dist_type const& dist) -> void {
		std::string sub_table = read_file(tex_subtable);
		pos = sub_table.find("^KR^");
		sub_table.replace(pos, 4, get_mu_str(mu));
		pos = sub_table.find("^COLS^");
		while (pos != std::string::npos) {
			sub_table.replace(pos, 6, ncols);
			pos = sub_table.find("^COLS^", pos);
		}
		std::string amps{};
		for (uint i=0; i<cols.size(); ++i)
			amps += " &";
		pos = sub_table.find("^AMPS^");
		sub_table.replace(pos, 6, amps);
	
		table_text += sub_table;
	
		for (uint i=0; i<dist.at(0).size(); ++i) {
			double x = xtab.at(i);
			table_text += scientificToLatex(x, 1, benchmark_format) + " & ";
				
			for (uint j=0; j<dist.size()-1; ++j)
				table_text += format_val(dist.at(j).at(i)) + " & ";
			table_text += format_val(dist.back().at(i));
			
			table_text += " \\\\\n";
		}
	};

	std::array<double, 3> mus{1.0, 0.5, 2.0};
	for (uint i=0; i<3; ++i)
		do_subtable(mus[i], arrs[i]);

	std::ifstream table_footer_s(tex_table_footer);
	std::string table_footer{std::istreambuf_iterator<char>(table_footer_s), std::istreambuf_iterator<char>{}};
	table_text += table_footer;

	fs::path latex_build_dir = fs::current_path()/"latex";
	if (!fs::exists(latex_build_dir)) {
		if (!fs::create_directory(latex_build_dir))
			log(LOG_ERROR, "util.hpp", "Failed to create latex build directory.");
		log(LOG_INFO, "util.hpp", "'latex' directory created.");
	} else {
	    log(LOG_INFO, "util.hpp", "'latex' directory exists. Continuing.");
	}

	std::string title = filename + ".tex";
	fs::path latex_file_path = latex_build_dir/title;
	std::ofstream latex_file(latex_file_path);
	latex_file << table_text;
	latex_file.close();

	std::string command = "pdflatex -interaction=batchmode -output-directory latex " + title;
	system(command.c_str());
	log(LOG_INFO, "util.hpp", "Cleaning up auxilliary files...");

	auto output_dir = (format == 0 || format == 2) ? "diffs" : "tables";
	fs::path output_dir_path = fs::current_path()/fs::path(output_dir);
	if (!fs::exists(output_dir_path)) {
		if (!fs::create_directory(output_dir_path))
			log(LOG_ERROR, "util.hpp", "Failed to create output directory '{}'.", output_dir_path.string());
	}
	fs::path pdf_path(fs::current_path()/fs::path("latex")/fs::path(filename + ".pdf"));
	log(LOG_INFO, "util.hpp", "{}  ->  {}", pdf_path.string(), output_dir_path.string());
	fs::copy(pdf_path, output_dir_path, fs::copy_options::overwrite_existing);
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
