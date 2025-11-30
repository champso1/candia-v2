#include "Candia-v2/Candia.hpp"
#include "Candia-v2/Common.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <filesystem>
#include <functional>

static char const* const TABLE_HEADER =
	"\\documentclass{article}\n"
	"\\usepackage{multicol}\n"
	"\\usepackage{graphicx} % Required for inserting images\n"
	"\\usepackage{a4}\n"
	"\\usepackage{amsmath}\n"
	"\\usepackage{amssymb}\n"
	"\\usepackage[\n"
	"  colorlinks=true\n"
	"  ,urlcolor=blue\n"
	"  ,anchorcolor=blue\n"
	"  ,citecolor=blue\n"
	"  ,filecolor=blue\n"
	"  ,linkcolor=blue\n"
	"  ,menucolor=blue\n"
	"  ,linktocpage=true\n"
	"  ,pdfproducer=medialab\n"
	"  ,pdfa=true\n"
	"]{hyperref}\n"
	"\\usepackage[capitalize]{cleveref}\n"
	"\\newcommand{\\TC}[1]{\\textcolor{ForestGreen}{\\bf[NOTE: TC -- #1]}}\n"
	"\\textheight 23.0cm \\textwidth 16.5cm\n"
	"\\oddsidemargin -0.1cm \\evensidemargin -0.1cm\n"
	"\\topmargin -1.5cm  % for hep-ph\n"
	"\\begin{document}\n"
	"\\begin{table}[htp]\n"
	"    \\caption{\n"
	"        Percentage error between Candia's results and MSHT's results,\n"
	"        both with the FHMRUVV parameterization\n"
	"    }\n"
	"    \\label{tab:n3lo_vfns_fhmruvv_msht}\n"
	"    \\begin{center}\n"
	"    \\vspace{5mm}\n"
	"    \\begin{tabular}{||c||r|r|r|r|r|r|r|r||}\n"
	"    \\hline \\hline\n"
	"    \\multicolumn{9}{||c||}{} \\\\[-3mm]\n"
	"    \\multicolumn{9}{||c||}{$\\, n_f = 3\\ldots 5\\,$,\n"
	"        $\\,\\mu_{\\rm f}^2 = 10^4 \\mbox{ GeV}^2$} \\\\\n"
	"    \\multicolumn{9}{||c||}{} \\\\[-0.3cm]\n"
	"    \\hline \\hline\n"
	"    \\multicolumn{9}{||c||}{} \\\\[-3mm]\n"
	"    \\multicolumn{1}{||c||}{$x$} &\n"
	"    \\multicolumn{1}{c|} {$xu_v$} &\n"
	"    \\multicolumn{1}{c|} {$xd_v$} &\n"
	"    \\multicolumn{1}{c|} {$xL_-$} &\n"
	"    \\multicolumn{1}{c|} {$xL_+$} &\n"
	"    \\multicolumn{1}{c|} {$xs_+$} &\n"
	"    \\multicolumn{1}{c|} {$xc_+$} &\n"
	"    \\multicolumn{1}{c|} {$xb_+$} &\n"
	"    \\multicolumn{1}{c||}{$xg$} \\\\[0.5mm]\n"
	"\\hline \\hline\n"
	"\\multicolumn{9}{||c||}{} \\\\[-3mm]\n"
	"\\multicolumn{9}{||c||}{$\\mu_{\\rm r}^2 = \\ \\mu_{\\rm f}^2$} \\\\\n"
	"\\multicolumn{9}{||c||}{} \\\\[-0.3cm]\n"
	"\\hline \\hline\n"
	" & & & & & & & \\\\[-0.3cm]\n";

static char const* const TABLE_FOOTER =
	"\\hline \\hline\n"
	"\\end{tabular}\n"
	"\\end{center}\n"
	"\\end{table}\n"
	"\\end{document}";



static std::string scientificToLatex(double num, uint precision)
{
	std::ostringstream ss{};
	ss << std::scientific << std::setprecision(precision) << num;
	std::string str = ss.str();

	auto e_pos = str.find("e");
	
	std::string mantissa = str.substr(0, e_pos);
	int exponent = std::stoi(str.substr(e_pos+1, std::string::npos));
	// ^ we convert to number then back to string here to remove leading zero

	std::ostringstream out_ss{};
	out_ss << mantissa << "$^{" << std::showpos << exponent << "}$";
	return out_ss.str();
}


static std::string percentToLatex(double num, uint precision)
{
	(void)precision;
	std::ostringstream ss{};
	ss << std::fixed << std::setprecision(2);
	if (num > 10.0)
		ss << std::scientific;
	ss << (num * 100.0) << "\\%";
	return ss.str();
}



namespace Candia2
{
	void DGLAPSolver::OutputLatexTable(MultiDimVector<double, 2>::type const& dists, std::string const& title, uint format)
	{
		namespace fs = std::filesystem;
		
		fs::path latex_build_dir{fs::current_path()};
		latex_build_dir /= "latex";
		if (!fs::exists(latex_build_dir))
		{
			if (!fs::create_directory(latex_build_dir))
			{
				std::cerr << "[ERROR] DGLAPSolver::OutputLatexTable(): failed to create latex build directory\n";
				exit(EXIT_FAILURE);
			}
		}

		fs::path latex_file_path{latex_build_dir};
		latex_file_path /= title;
		
		std::ofstream latex_file(latex_file_path);
		if (!latex_file)
		{
			std::cerr << "[ERROR] DGLAPSolver::OutputLatexTable(): failed to open "
					  << std::quoted(latex_file_path.string())
					  << "\n";
			exit(EXIT_FAILURE);
		}

		latex_file << TABLE_HEADER;

		
		// conventions:
		// q_v = q - qbar
		// L- = (dbar - ubar)
		// L+ = 2(dbar + ubar)
		// q+ - q + qbar
		//
		// 0      gluons         g
		// 1-6    quarks         u,d,s,c,b,t
		// 7-12   antiquarks     au,ad,as,ac,ab,at
		//
		// the list contains:
		// xu_v  xd_v  xL-  xL+  xs_v  xs_+  xc_+  xg

		std::function<std::string(double,uint)> format_func;
		switch (format)
		{
			case PERCENT: format_func = percentToLatex; break;
			case SCIENTIFIC: format_func = scientificToLatex; break;
			default:
			{
				std::cerr << "[ERROR] DGLAPSolver::OutputLatexTable(): invalid format "
						  << std::quoted(std::to_string(format))
						  << '\n';
				exit(EXIT_FAILURE);
			}
		}

		std::vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9};
		for (uint ix=0; ix<xtab.size(); ++ix)
		{
			double x = xtab.at(ix);
			latex_file << scientificToLatex(x, 1) << " & "; // always scientific
				
			for (uint i=0; i<dists.size()-1; ++i)
				latex_file << format_func(dists.at(i).at(ix), 4) << " & ";
			latex_file << format_func(dists.at(7).at(ix), 4);
			
			latex_file << " \\\\\n";
		}
		
		latex_file << TABLE_FOOTER;
		latex_file.close();

		std::string command = "pdflatex -output-directory latex " + title;
		system(command.c_str());
	}
} // namespace Candia2



// old function just in case
/*
void DGLAPSolver::OutputLatexTable(MultiDimVector<double, 2>::type const& dists)
{
	namespace fs = std::filesystem;
	fs::path latex_build_dir{fs::current_path()};
	latex_build_dir /= "latex";
	if (!fs::exists(latex_build_dir))
	{
		if (!fs::create_directory(latex_build_dir))
		{
			std::cerr << "[ERROR] DGLAPSolver::OutputLatexTable(): failed to create latex build directory\n";
			exit(EXIT_FAILURE);
		}
	}

	fs::path latex_file_path{latex_build_dir};
	latex_file_path /= "table.tex";
		
	std::ofstream latex_file(latex_file_path);
	if (!latex_file)
	{
		std::cerr << "[ERROR] DGLAPSolver::OutputLatexTable(): failed to open " << std::quoted("table.tex") << "\n";
		exit(EXIT_FAILURE);
	}

	latex_file << TABLE_HEADER;

	std::vector<double>::size_type size = dists.at(0).size();
		
	// conventions:
	// q_v = q - qbar
	// L- = (dbar - ubar)
	// L+ = 2(dbar + ubar)
	// q+ - q + qbar
	//
	// 0      gluons         g
	// 1-6    quarks         u,d,s,c,b,t
	// 7-12   antiquarks     au,ad,as,ac,ab,at
	//
	// the list contains:
	// xu_v  xd_v  xL-  xL+  xs_v  xs_+  xc_+  xg
	std::array<std::vector<double>, 8> table_dists = [s = size]{
		std::array<std::vector<double>, 8> temp;
		for (auto& v : temp)
			v = std::vector<double>(s);
		return temp;
	}();

		
	for (uint k=0; k<size; ++k)
	{
		table_dists.at(0).at(k) = dists[1][k] - dists[1+6][k];
		table_dists.at(1).at(k) = dists[2][k] - dists[2+6][k];
		table_dists.at(2).at(k) = dists[2+6][k] - dists[1+6][k];
		table_dists.at(3).at(k) = 2.0*(dists[2+6][k] + dists[1+6][k]);
		table_dists.at(4).at(k) = dists[3][k] + dists[3+6][k];
		table_dists.at(5).at(k) = dists[4][k] + dists[4+6][k];
		table_dists.at(6).at(k) = dists[5][k] + dists[5+6][k];
		table_dists.at(7).at(k) = dists[0][k];
	}

	for (uint ix=0; ix<ntab.size()-2; ++ix)
	{
		double x = xtab.at(ix);
		latex_file << scientificToLatex(x, 1) << " & ";

		uint idx = ntab.at(ix);
				
		for (uint i=0; i<table_dists.size()-1; ++i)
			latex_file << scientificToLatex(table_dists.at(i).at(idx), 4) << " & ";
		latex_file << scientificToLatex(table_dists.at(7).at(idx), 4);
		latex_file << " \\\\\n";
	}
		
	latex_file << TABLE_FOOTER;
	latex_file.close();

	system("pdflatex -output-directory latex table.tex");

	// for convenience, copy the produced pdf file to the main output directory (i.e. out from latex/)
		
	fs::copy_file("build/table.pdf", "table.pdf");
}
*/
