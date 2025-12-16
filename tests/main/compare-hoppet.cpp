#include <cstdlib>
#include <ios>
#include <print>
#include <limits>
#include <sstream>
#include <string>
#include <algorithm>
#include <ranges>
#include <vector>
#include <filesystem>
#include <fstream>
using namespace std;
namespace fs = filesystem;

#include "Candia-v2/Common.hpp"
using namespace Candia2;

constexpr static char const* HOPPET_DATA_FILE{"out-nnlo.dat"};

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
	"        Percentage error between Candia's results and HOPPET's results.\n"
	"    }\n"
	"    \\begin{center}\n"
	"    \\vspace{5mm}\n"
	"    \\begin{tabular}{||c||r|r|r|r|r|r|r|r|r|r|r|}\n"
	"    \\hline \\hline\n"
	"    \\multicolumn{12}{||c||}{} \\\\[-3mm]\n"
	"    \\multicolumn{12}{||c||}{$\\, n_f = 3\\ldots 5\\,$,\n"
	"        $\\,\\mu_{\\rm f}^2 = 10^4 \\mbox{ GeV}^2$} \\\\\n"
	"    \\multicolumn{12}{||c||}{} \\\\[-0.3cm]\n"
	"    \\hline \\hline\n"
	"    \\multicolumn{12}{||c||}{} \\\\[-3mm]\n"
	"    \\multicolumn{1}{||c||}{$x$} &\n"
	"    \\multicolumn{1}{c|} {$xg$} &\n"
	"    \\multicolumn{1}{c|} {$xu$} &\n"
	"    \\multicolumn{1}{c|} {$xd$} &\n"
	"    \\multicolumn{1}{c|} {$xs$} &\n"
	"    \\multicolumn{1}{c|} {$xc$} &\n"
	"    \\multicolumn{1}{c|} {$xb$} &\n"
	"    \\multicolumn{1}{c|} {$xub$} &\n"
	"    \\multicolumn{1}{c|} {$xdb$} &\n"
	"    \\multicolumn{1}{c|} {$xsb$} &\n"
	"    \\multicolumn{1}{c|} {$xcb$} &\n"
	"    \\multicolumn{1}{c||}{$xbb$} \\\\[0.5mm]\n";

static char const* TABLE_SUBHEADER =
	"\\hline \\hline\n"
	"\\multicolumn{12}{||c||}{} \\\\[-3mm]\n"
	"\\multicolumn{12}{||c||}{$\\mu_{\\rm r}^2 = \\ %KR%\\mu_{\\rm f}^2$} \\\\\n"
	"\\multicolumn{12}{||c||}{} \\\\[-0.3cm]\n"
	"\\hline \\hline\n"
	" & & & & & & & & & & & \\\\[-0.3cm]\n";


static char const* const TABLE_FOOTER =
	"\\hline \\hline\n"
	"\\end{tabular}\n"
	"\\end{center}\n"
	"\\end{table}\n"
	"\\end{document}";

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
static std::string percentToLatex(double num, uint precision)
{
	(void)precision;
	ostringstream ss{};
	ss << fixed << setprecision(2);
	if (num > 10.0)
		ss << scientific;
	ss << (num * 100.0) << "\\%";
	return ss.str();
}

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


void outputLatexTable(vector<vector<double>> const& diffs, vector<double> const& xtab, string const& basename)
{
	fs::path latex_build_dir{fs::current_path()/"latex"};
	if (!fs::exists(latex_build_dir))
	{
		if (!fs::create_directory(latex_build_dir))
		{
			println("failed to create latex build directory");
			exit(EXIT_FAILURE);
		}
		println("\"latex\" directory created.");
	}
	else
		println("\"latex\" directory already exists.");
	
	string title = basename + ".tex";
	fs::path latex_file_path{latex_build_dir/title};
	ofstream latex_file(latex_file_path);

	latex_file << TABLE_HEADER;

	println("Table header has been written.");

	string subheader = tableSubheader(1.0);
	latex_file << subheader;
	for (uint ix=0; ix<xtab.size(); ++ix)
	{
		double x = xtab.at(ix);
		latex_file << scientificToLatex(x, 1) << " & ";
				
		for (uint i=0; i<diffs.size()-1; ++i)
			latex_file << percentToLatex(diffs.at(i).at(ix), 4) << " & ";
		latex_file << percentToLatex(diffs.back().at(ix), 4);
			
		latex_file << " \\\\\n";
	}
	println("Done.");	
		
	latex_file << TABLE_FOOTER;
	latex_file.close();

	print("Wrote table header. Running pdflatex... ");
	string command = "pdflatex -interaction=batchmode -output-directory latex " + title;
	system(command.c_str());
	println("Cleaning up auxilliary files...");

	fs::path pdf_path{fs::current_path()/"latex"}, new_pdf_path{fs::current_path()};
	ranges::filter_view tex_pdf_files =
		fs::directory_iterator{fs::current_path()/"latex"}
		| ranges::views::filter([&b = basename](fs::directory_entry const& e)
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
		println("[ERROR] Must provide a file to read from.");
		exit(EXIT_FAILURE);
	}

	fs::path candia_data_path{argv[1]};
	if (!fs::exists(candia_data_path))
	{
	    println("[ERROR] Failed to find file \"{}\"", candia_data_path.string());
		exit(EXIT_FAILURE);
	}

	fs::path hoppet_data_path{HOPPET_DATA_FILE};
	if (!fs::exists(hoppet_data_path))
	{
		println("[ERROR] Failed to find existing HOPPET data in file \"{}\"", hoppet_data_path.string());
		exit(EXIT_FAILURE);
	}

	// fill vectors with zeros
	// 11 dists (all of them minus top quark)
	// and 9 tabulated x values each (1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.9)
	MultiDimVector<double, 2>::type candia_data(11, std::vector<double>(9, 0.0));
	MultiDimVector<double, 2>::type	hoppet_data{candia_data};
	
	// output vector for relative diffs
	MultiDimVector<double, 2>::type diffs{candia_data};
	
	// read from the candia file
	ifstream candia_stream{candia_data_path};
	candia_stream.ignore(numeric_limits<streamsize>::max(), '\n'); // ignore the comment line

	// contains xtab and ntab lines, read those in
	vector<double> xtab{};
	vector<int> ntab{};
	double temp1{};
	int temp2;
	string line{};
	
	getline(candia_stream, line);
	istringstream iss{line};
	while (iss >> temp1)
		xtab.push_back(temp1);
	
	getline(candia_stream, line);
	iss = istringstream{line};
	while (iss >> temp2)
		ntab.push_back(temp2);

	// grab everything into a table F, we will parse this to match after
	MultiDimVector<double, 2>::type F(13, std::vector<double>{});
	double x{};
	while (getline(candia_stream, line))
	{
		iss = istringstream{line};
		iss >> x; // we ignore this since we have the tabulated x values already
		for (int i=0; i<F.size(); ++i)
		{
			iss >> temp1;
			F.at(i).push_back(temp1);
		}
	}

	// we skip the 1.0 from the end of the candia-v2 distributions
	for (uint ik=0; ik<ntab.size()-1; ++ik)
	{
		int k = ntab.at(ik);
		candia_data.at(0).at(ik) =  F[0][k];
		candia_data.at(1).at(ik) =  F[1][k];
		candia_data.at(2).at(ik) =  F[2][k];
		candia_data.at(3).at(ik) =  F[3][k];
		candia_data.at(4).at(ik) =  F[4][k];
		candia_data.at(5).at(ik) =  F[5][k];
		candia_data.at(6).at(ik) =  F[6+1][k];
		candia_data.at(7).at(ik) =  F[6+2][k];
		candia_data.at(8).at(ik) =  F[6+3][k];
		candia_data.at(9).at(ik) =  F[6+4][k];
		candia_data.at(10).at(ik) = F[6+5][k];
	}

	ifstream hoppet_stream{hoppet_data_path};
	if (!hoppet_stream)
	{
	    println("[ERROR] Failed to open a file stream for the HOPPET data file \"{}\"", HOPPET_DATA_FILE);
		exit(EXIT_FAILURE);
	}
	hoppet_stream.ignore(numeric_limits<streamsize>::max(), '\n'); // comment line

	uint idx = 0;
	double temp{};
	while (getline(hoppet_stream, line))
	{
		iss = istringstream{line};
		iss >> temp; // don't actually care about the x value, we already know them
		for (int j=0; j<hoppet_data.size(); ++j)
			iss >> hoppet_data.at(j).at(idx);
		idx++;
	}

	if ((hoppet_data.size() != candia_data.size()) || (hoppet_data.at(0).size() != candia_data.at(0).size()))
	{
		println("Hoppet and Candia data size mismatch");
		exit(EXIT_FAILURE);
	}

	auto reldiff =
		[](double candia, double hoppet) -> double {
			return abs(candia-hoppet)/hoppet;
		};

	for (uint j=0; j<candia_data.size(); ++j)
	{
		for (uint k=0; k<candia_data.at(0).size(); ++k)
		{
			double candia = candia_data.at(j).at(k);
			double hoppet = hoppet_data.at(j).at(k);
			diffs.at(j).at(k) = reldiff(candia, hoppet);
		}
	}
	
    ostringstream filename_ss{};
	filename_ss << "diffs-hoppet-" << candia_data_path.filename().string().substr(0, candia_data_path.filename().string().rfind('.'));
	println("Outputting to filename \"{}\"", filename_ss.str());
	outputLatexTable(diffs, vector<double>(xtab.begin(), xtab.end()-1), filename_ss.str());
}
