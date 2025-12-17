#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <filesystem>
#include <fstream>
#include <print>
#include <unordered_map>
using namespace std;
using uint = unsigned;
namespace fs = filesystem;
using DistVec = vector<vector<double>>;



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
	"    \\multicolumn{1}{c||}{$xg$} \\\\[0.5mm]\n";

static char const* const TABLE_HEADER_CANDIAV1 =
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
	"    \\multicolumn{3}{||c||}{} \\\\[-3mm]\n"
	"    \\multicolumn{3}{||c||}{$\\, n_f = 3\\ldots 5\\,$,\n"
	"        $\\,\\mu_{\\rm f}^2 = 10^4 \\mbox{ GeV}^2$} \\\\\n"
	"    \\multicolumn{3}{||c||}{} \\\\[-0.3cm]\n"
	"    \\hline \\hline\n"
	"    \\multicolumn{3}{||c||}{} \\\\[-3mm]\n"
	"    \\multicolumn{1}{||c||}{$x$} &\n"
	"    \\multicolumn{1}{c|} {$xg$} &\n"
	"    \\multicolumn{1}{c||}{$xq^{(-)}$} \\\\[0.5mm]\n";

static char const* TABLE_SUBHEADER =
	"\\hline \\hline\n"
	"\\multicolumn{9}{||c||}{} \\\\[-3mm]\n"
	"\\multicolumn{9}{||c||}{$\\mu_{\\rm r}^2 = \\ %KR%\\mu_{\\rm f}^2$} \\\\\n"
	"\\multicolumn{9}{||c||}{} \\\\[-0.3cm]\n"
	"\\hline \\hline\n"
	" & & & & & & & \\\\[-0.3cm]\n";

static char const* TABLE_SUBHEADER_CANDIAV1 =
	"\\hline \\hline\n"
	"\\multicolumn{3}{||c||}{} \\\\[-3mm]\n"
	"\\multicolumn{3}{||c||}{$\\mu_{\\rm r}^2 = \\ %KR%\\mu_{\\rm f}^2$} \\\\\n"
	"\\multicolumn{3}{||c||}{} \\\\[-0.3cm]\n"
	"\\hline \\hline\n"
	" & & \\\\[-0.3cm]\n";


static char const* const TABLE_FOOTER =
	"\\hline \\hline\n"
	"\\end{tabular}\n"
	"\\end{center}\n"
	"\\end{table}\n"
	"\\end{document}";


static string tableSubheader(double _kr)
{
	string kr{};
	if (_kr != 1.0)
	{
		ostringstream oss{};
		oss << setprecision(1) << _kr;
		kr = oss.str();
	}

	string _subheader{TABLE_SUBHEADER};
    string find_expr{"%KR%"};
	string subheader = _subheader.replace(_subheader.find(find_expr), find_expr.size(), kr);
	return subheader;
}





static string scientificToLatex(double num, uint precision)
{
	ostringstream ss{};
	ss << scientific << setprecision(precision) << num;
	string str = ss.str();

	auto e_pos = str.find("e");
	
	string mantissa = str.substr(0, e_pos);
	int exponent = stoi(str.substr(e_pos+1, string::npos));
	// ^ we convert to number then back to string here to remove leading zero

	ostringstream out_ss{};
	out_ss << mantissa << "$^{" << showpos << exponent << "}$";
	return out_ss.str();
}


std::tuple<DistVec, vector<int>, vector<double>> readDatafile(fs::path path)
{
	print("Reading data from file \"{}\"... ", path.filename().string());
	
	ifstream file_stream{path};
	file_stream.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore the comment line

	vector<double> xtab{};
	vector<int> ntab{};
	
	double temp1{};
	int temp2{};
	string line{};

	// read in xtab array
	getline(file_stream, line);
	istringstream iss{line};
	while (iss >> temp1)
		xtab.push_back(temp1);

	// read in ntab array
	getline(file_stream, line);
	iss = istringstream{line};
	while (iss >> temp2)
		ntab.push_back(temp2);

	// read in rest of data points
	vector<double> X{};
	DistVec F{};
	F.resize(13);
	while (getline(file_stream, line))
	{
		iss = istringstream{line};
		iss >> temp1;
		X.push_back(temp1);
		for (int i=0; i<F.size(); ++i)
		{
			iss >> temp1;
			F.at(i).push_back(temp1);
		}
	}

	println("Done.");
	println("Assigning new distributions in accordance with the table... ");

	// must create required dists from the regular ones
	DistVec dists(8, vector<double>(11, 0.0));
	double size = ntab.size();
	for (uint ik=0; ik<size-1; ++ik) // -1 because we don't want x=1.0
	{
		int k = ntab.at(ik);
		
		dists.at(0).at(ik) = F[1][k] - F[1+6][k];
		dists.at(1).at(ik) = F[2][k] - F[2+6][k];
		dists.at(2).at(ik) = F[2+6][k] - F[1+6][k];
		dists.at(3).at(ik) = 2.0*(F[2+6][k] + F[1+6][k]);
		dists.at(4).at(ik) = F[3][k] + F[3+6][k];
		dists.at(5).at(ik) = F[4][k] + F[4+6][k];
		dists.at(6).at(ik) = F[5][k] + F[5+6][k];
		dists.at(7).at(ik) = F[0][k];

		println("x={}:", xtab.at(ik));
		println("    xv: {} - {}", F[1][k], F[1+6][k]);
		println("    dv: {} - {}", F[2][k], F[2+6][k]);
		println("    L-: {} - {}", F[2+6][k], F[1+6][k]);
		println("    L+: 2({} + {})", F[2+6][k], F[1+6][k]);
		println("    s+: {} + {}", F[3][k], F[3+6][k]);
		println("    c+: {} + {}", F[4][k], F[4+6][k]);
		println("    b+: {} + {}", F[5][k], F[5+6][k]);
		println("    g:  {}", F[0][k]);
	}
	println("Done.");
	return {dists, ntab, xtab};
}

std::tuple<DistVec, vector<int>, vector<double>> readDatafile_Candiav1(fs::path path)
{
	print("Reading data from file \"{}\"... ", path.filename().string());
	
	ifstream file_stream{path};
	file_stream.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore the comment line

	vector<double> xtab{};
	vector<int> ntab{};
	
	double temp1{};
	int temp2{};
	string line{};

	// read in xtab array
	getline(file_stream, line);
	istringstream iss{line};
	while (iss >> temp1)
		xtab.push_back(temp1);

	// read in ntab array
	getline(file_stream, line);
	iss = istringstream{line};
	while (iss >> temp2)
		ntab.push_back(temp2);

	// read in rest of data points
	vector<double> X{};
	DistVec F{};
	F.resize(13);
	while (getline(file_stream, line))
	{
		iss = istringstream{line};
		iss >> temp1;
		X.push_back(temp1);
		for (int i=0; i<F.size(); ++i)
		{
			iss >> temp1;
			F.at(i).push_back(temp1);
		}
	}

	println("Done.");
	println("Assigning new distributions in accordance with the table... ");

	// must create required dists from the regular ones
	DistVec dists(2, vector<double>(9, 0.0));
	double size = ntab.size();
	for (uint ik=0; ik<size-1; ++ik) // -1 because we don't want x=1.0
	{
		int k = ntab.at(ik);
		
		dists.at(0).at(ik) = F[0][k];
		double res = 0.0;
		for (uint j=1; j<=6; ++j)
			res += F[j][k] - F[j+6][k];
		dists.at(1).at(ik) = res;
	}
	println("Done.");
	return {dists, ntab, xtab};
}

void outputLatexTable(unordered_map<double, fs::path> const& paths, string const& basename)
{
	fs::path latex_build_dir{fs::current_path()/"latex"};
	if (!fs::exists(latex_build_dir))
	{
		if (!fs::create_directory(latex_build_dir))
		{
			cerr << "[ERROR] DGLAPSolver::OutputLatexTable(): failed to create latex build directory\n";
			exit(EXIT_FAILURE);
		}
		println("\"latex\" directory created.");
	}
	else
		println("\"latex\" directory already exists.");
	
	string title = basename + ".tex";
	fs::path latex_file_path{latex_build_dir/title};
	ofstream latex_file(latex_file_path);
	if (!latex_file)
	{
		cerr << "[ERROR] read.cpp: failed to open "
				  << quoted(latex_file_path.string())
				  << "\n";
		exit(EXIT_FAILURE);
	}

	latex_file << TABLE_HEADER;

	println("Table header has been written.");

		
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

    for (pair<double, fs::path> const& x : paths)
	{
		auto [dists, ntab, xtab] = readDatafile(x.second);
		print("Writing data from file \"{}\"... ", x.second.filename().string());
		string subheader = tableSubheader(x.first);
		latex_file << subheader;
		for (uint ix=0; ix<xtab.size()-1; ++ix)
		{
			double x = xtab.at(ix);
			latex_file << scientificToLatex(x, 1) << " & ";
				
			for (uint i=0; i<dists.size()-1; ++i)
				latex_file << scientificToLatex(dists.at(i).at(ix), 4) << " & ";
			latex_file << scientificToLatex(dists.at(7).at(ix), 4);
			
			latex_file << " \\\\\n";
		}
		println("Done.");
	}
	
	
		
	latex_file << TABLE_FOOTER;
	latex_file.close();

	print("Wrote table header. Running pdflatex... ");
	string command = "pdflatex -interaction=batchmode -output-directory latex " + title;
	system(command.c_str());
	println("Cleaning up auxilliary files...");

	fs::path pdf_path{fs::current_path()/"latex"}, new_pdf_path{fs::current_path()};
	ranges::filter_view tex_pdf_files =
		fs::directory_iterator{fs::current_path()/"latex"}
		| ranges::views::filter([b = basename](fs::directory_entry const& e)
		{
			string filename = e.path().filename().string();
		    return filename.contains(b);
		});

	ranges::for_each(tex_pdf_files, [](fs::directory_entry const& e) {
		fs::path path = e.path();
	    string filename = path.filename().string();
		string ext = filename.substr(filename.rfind('.')+1, string::npos);
		if (ext.compare("pdf") == 0)
		{
			fs::path new_pdf_path{fs::current_path()/filename};	
			fs::copy(path, new_pdf_path, fs::copy_options::overwrite_existing);
			fs::remove(path);
		}
		else if (ext.compare("tex") != 0)
			fs::remove(path);
	});
}


int main(int argc, char *argv[])
{
	if (argc != 2)
	{
	    println("[ERROR] read.cpp: Must provide an evolution to read from.");
		exit(EXIT_FAILURE);
	}

	string basename{argv[1]};

	// find all of the associated files
	unordered_map<double, fs::path> datafiles{};
	for (fs::directory_entry const& entry : fs::directory_iterator{fs::current_path()})
	{
		string _path = entry.path().filename().string();
		string _ext = _path.substr(_path.rfind('.')+1, string::npos);
		if ((_path.find(basename) == 0) && (_ext.compare("dat") == 0))
		{
			streamsize r_pos = _path.find("-r");
			streamsize dot_pos = _path.rfind('.');
			if (r_pos == string::npos || dot_pos == string::npos)
			{
				println("[ERROR] read.cpp: Failed to parse the evolution given by: \"{}\".", _path);
				exit(EXIT_FAILURE);
			}
			string kr_extract = _path.substr(r_pos + 2, dot_pos-r_pos-2);
		    double kr = stold(kr_extract);
			datafiles[kr] = entry.path();
		}
	}
	println("Found {} datafiles in the current directory corresponding to the evolution \"{}\":", datafiles.size(), basename);
	ranges::for_each(datafiles, [](std::pair<double, fs::path> const& x) -> void {
		println("    {:.1} -> \"{}\".", x.first, x.second.filename().string());
	});
	if (datafiles.size() == 0)
	{
		println(cerr, "No datafiles found. Exiting...");
		exit(EXIT_FAILURE);
	}

	outputLatexTable(datafiles, basename);
}
