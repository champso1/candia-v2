#include "Candia-v2/Grid.hpp"
#include "Candia-v2/SplittingFn.hpp"

using SplitFunc = Candia2::SplittingFunction;

#include <vector>
#include <filesystem>
#include <fstream>
#include <iomanip>
using namespace std;

int main()
{
	vector<double> xtab{1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Candia2::Grid grid(xtab, 401);

	Candia2::SplittingFunction::Update(4, 8.3333333);

	Candia2::P1qq p0ps{};
	Candia2::P1nsm p1ps{};
	Candia2::P2ps p2ps{};
	Candia2::P3ps p3ps{};
	Candia2::P3ps p3ps_a1(1);
	Candia2::P3ps p3ps_a2(2);
	Candia2::P2nsm p2nsm{};
	Candia2::P2nsp p2nsp{};
	Candia2::P3nsm p3nsm{};
	Candia2::P3nsm p3nsm_a1(1);
	Candia2::P3nsm p3nsm_a2(2);
	Candia2::P3nsp p3nsp{};
	Candia2::P3nsp p3nsp_a1(1);
	Candia2::P3nsp p3nsp_a2(2);

	std::vector<double> x_vals{};
	std::vector<double> p0ps_vals{}, p1ps_vals{}, p2ps_vals{}, p3ps_vals{};
	std::vector<double> p3ps_a1_vals{}, p3ps_a2_vals{};
	std::vector<double> p2nsm_vals{}, p2nsp_vals{};
	std::vector<double> p3nsm_vals{}, p3nsp_vals{};
	std::vector<double> p3nsm_a1_vals{}, p3nsm_a2_vals{}, p3nsp_a1_vals{}, p3nsp_a2_vals{};

	for (uint k=0; k<=grid.Size()-1; k++)
	{
		double x = grid.At(k);
		x_vals.push_back(x);

		p0ps_vals.push_back(x*p0ps.Regular(x) * 2.0);
		p1ps_vals.push_back(x*p1ps.Regular(x) * 4.0);
		p2ps_vals.push_back(x*p2ps.Regular(x)/2000.0 * 8.0);
		
		p3ps_vals.push_back(x*p3ps.Regular(x)/25000.0 * 16.0);
		p3ps_a1_vals.push_back(x*p3ps_a1.Regular(x)/25000.0 * 16.0);
		p3ps_a2_vals.push_back(x*p3ps_a2.Regular(x)/25000.0 * 16.0);

		p2nsm_vals.push_back(8.0*(1.0-x)*(p2nsm.Regular(x) + p2nsm.Plus(x)/(1.0-x))/2000.0);
		p2nsp_vals.push_back(8.0*(1.0-x)*(p2nsp.Regular(x) + p2nsp.Plus(x)/(1.0-x))/2000.0);
		
		p3nsm_vals.push_back(16.0*(1.0-x)*(p3nsm.Regular(x) + p3nsm.Plus(x)/(1.0-x))/25000.0);
		p3nsp_vals.push_back(16.0*(1.0-x)*(p3nsp.Regular(x) + p3nsp.Plus(x)/(1.0-x))/25000.0);
		
		p3nsm_a1_vals.push_back(16.0*(1.0-x)*(p3nsm_a1.Regular(x) + p3nsm_a1.Plus(x)/(1.0-x))/25000.0);
		p3nsm_a2_vals.push_back(16.0*(1.0-x)*(p3nsm_a2.Regular(x) + p3nsm_a2.Plus(x)/(1.0-x))/25000.0);
		p3nsp_a1_vals.push_back(16.0*(1.0-x)*(p3nsp_a1.Regular(x) + p3nsp_a1.Plus(x)/(1.0-x))/25000.0);
		p3nsp_a2_vals.push_back(16.0*(1.0-x)*(p3nsp_a2.Regular(x) + p3nsp_a2.Plus(x)/(1.0-x))/25000.0);
	}

	filesystem::path outfile_path{"ps.dat"};
	ofstream outfile{outfile_path};
	if (!outfile)
		exit(5);

    for (uint k=0; k<=grid.Size()-1; k++)
	{
	    outfile << scientific << setprecision(9);
		outfile << x_vals[k] << '\t';

		outfile << fixed << setprecision(9);
		outfile << p1ps_vals[k] << '\t';
		outfile << p2ps_vals[k] << '\t';
		outfile << p3ps_vals[k] << '\t';
		outfile << p3ps_a1_vals[k] << '\t';
		outfile << p3ps_a2_vals[k] << '\t';
	}
	outfile << endl;

	outfile.close();

    outfile_path = "ns.dat";
	outfile.open(outfile_path);
	if (!outfile)
		exit(5);

    for (uint k=0; k<=grid.Size()-1; k++)
	{
	    outfile << scientific << setprecision(9);
		outfile << x_vals[k] << '\t';

		outfile << fixed << setprecision(9);
		outfile << p2nsm_vals[k] << '\t';
		outfile << p2nsp_vals[k] << '\t';
		outfile << p3nsm_vals[k] << '\t';
		outfile << p3nsp_vals[k] << '\t';
		outfile << p3nsm_a1_vals[k] << '\t';
		outfile << p3nsm_a2_vals[k] << '\t';
		outfile << p3nsp_a1_vals[k] << '\t';
		outfile << p3nsp_a2_vals[k] << '\n';
	}
	outfile << endl;
	
	outfile.close();
}
