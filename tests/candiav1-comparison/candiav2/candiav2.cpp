#include <iomanip>
#include <iostream>
#include <string>
#include <limits>
using namespace std;

#include "Candia-v2/Grid.hpp"
#include "Candia-v2/AlphaS.hpp"
#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/Distribution.hpp"
#include "Candia-v2/Candia.hpp"
using namespace Candia2;


const static string output_prefix = "output/output-candiav2/";
using SplitFunc = shared_ptr<SplittingFunction>;
void OutputFileSplitFunc(string const& filename, Grid const& grid, vector<SplitFunc> const& funcs);
void OutputFileAlphaSBetas(string const& filepath, AlphaS &as);
void OutputFileInitialDistributions(string const& filepath, Grid const& grid, shared_ptr<Distribution> dist);
void OutputFileGauleg(string const& filepath, Grid const& grid);
void OutputFileConvolution(string const& filepath, vector<double> const& vec);
void OutputFileInterpolation(string const& filepath, vector<double> const& vec);


int main() {
    // define the grid points
	vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
	Grid grid(xtab, 801);

	shared_ptr<Distribution> dist = make_shared<LesHouchesDistribution>();
	AlphaS alpha_s(0, dist->Q0(), dist->Alpha0(), dist->Masses()); // order here doesn't matter, we only need coeffs
	alpha_s.CalculateThresholdValues(100.0);
	alpha_s.Update(4); // simply set to four flavors
	SplittingFunction::Update(4, alpha_s.Beta0());

	vector<SplitFunc> LO
	{
		make_shared<P0ns>(),
		make_shared<P0qq>(),
		make_shared<P0qg>(),
		make_shared<P0gq>(),
		make_shared<P0gg>()
	};

	vector<SplitFunc> NLO
	{
		make_shared<P1nsp>(),
		make_shared<P1nsm>(),
		make_shared<P1qq>(),
		make_shared<P1qg>(),
		make_shared<P1gq>(),
		make_shared<P1gg>()
	};

	vector<SplitFunc> NNLO
	{
		make_shared<P2nsp>(),
		make_shared<P2nsm>(),
		make_shared<P2nsv>(),
		make_shared<P2qq>(),
		make_shared<P2qg>(),
		make_shared<P2gq>(),
		make_shared<P2gg>()
	};

	vector<SplitFunc> N3LO
	{
		make_shared<P3nsp>(),
		make_shared<P3nsm>(),
		make_shared<P3nsv>()
	};


	OutputFileSplitFunc("out-LO-splitfuncs.dat", grid, LO);
	OutputFileSplitFunc("out-NLO-splitfuncs.dat", grid, NLO);
	OutputFileSplitFunc("out-NNLO-splitfuncs.dat", grid, NNLO);
	OutputFileSplitFunc("out-N3LO-splitfuncs.dat", grid, N3LO);

	OutputFileAlphaSBetas("out-alphas_betas.dat", alpha_s);
	OutputFileInitialDistributions("out-init_dists.dat", grid, dist);
	OutputFileGauleg("out-gauleg.dat", grid);

	// this test is to fill a vector with 1's
	// and convolute it with a non-trivial splitting function
	// to test if the convolution is correct
	vector<double> vec(801, 1.0);
	vector<double> out(801);
	std::shared_ptr<SplittingFunction> p1nsp = make_shared<P1nsp>();

	for (uint k=0; k<grid.Size(); k++)
	{
		out[k] = grid.Convolution(vec, p1nsp, k);
	}
	OutputFileConvolution("out-conv.dat", out);

	
	vector<double> pts(801);
	out.clear();
	out.resize(67);
	double x;
	for (uint i=0; i<801; i++)
	{
		x = grid[i];
		pts[i] = sin(x);
	}

	for (uint i=0; i<67; i++)
	{
		x = static_cast<double>(i)/67.0;
		out[i] = grid.Interpolate(pts, x, true);
	}
    OutputFileInterpolation("out-interp.dat", out);



	DGLAPSolver solver(0, grid, 100.0, dist);
	solver.Evolve();
}






void OutputFileSplitFunc(string const& filename, Grid const& grid, vector<SplitFunc> const& funcs)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix.c_str(), filename.c_str());
	FILE * f = fopen(filename_buf, "w");

	SplitFunc func;
	double x;
	for (uint k=0; k<grid.Size(); k++)
	{
		x = grid[k];
		fprintf(f, "%15.9g\t", x);
		for (auto func : funcs)
		{
			fprintf(f, "%8.4lf\t", func->Regular(x));
			fprintf(f, "%8.4lf\t", func->Plus(x));
			fprintf(f, "%8.4lf\t", func->Delta(x));
		}
		fprintf(f, "\n");
	}
	
	fclose(f);
}

void OutputFileAlphaSBetas(string const& filepath, AlphaS &as)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix.c_str(), filepath.c_str());
	FILE * f = fopen(filename_buf, "w");

	for (uint nf=3; nf<=5; nf++)
	{
		as.Update(nf);
		double beta0 = as.Beta0();
		double beta1 = as.Beta1();
		double beta2 = as.Beta2();
		fprintf(f, "%d\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\n",
				nf, beta0, beta1, beta2, as.Pre(nf), as.Post(nf));
	}
	
	fclose(f);
}


void OutputFileInitialDistributions(string const& filepath, Grid const& grid, shared_ptr<Distribution> dist)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix.c_str(), filepath.c_str());
	FILE * f = fopen(filename_buf, "w");

	double x;
	for (uint k=0; k<grid.Size()-1; k++)
	{
		x = grid[k];
		fprintf(f, "%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\t%8.4lf\n",
				x, dist->xuv(x), dist->xdv(x),dist->xg(x), dist->xdb(x), dist->xub(x),dist-> xs(x));
	}
	
	fclose(f);
}

void OutputFileGauleg(string const& filepath, Grid const& grid)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix.c_str(), filepath.c_str());
	FILE * f = fopen(filename_buf, "w");

	
	
	for (uint i=0; i<30; i++)
	{
		fprintf(f, "%8.4lf\t%8.4lf\n", grid.Abscissae(i), grid.Weights(i));
	}
	
	fclose(f);
}

void OutputFileConvolution(string const& filepath, vector<double> const& vec)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix.c_str(), filepath.c_str());
	FILE * f = fopen(filename_buf, "w");

	for (const double x : vec)
	{
		fprintf(f, "%8.4lf\n", x);
	}
	
	fclose(f);
}

void OutputFileInterpolation(string const& filepath, vector<double> const& vec)
{
	char filename_buf[64] = {0};
	sprintf(filename_buf, "%s%s", output_prefix.c_str(), filepath.c_str());
	FILE * f = fopen(filename_buf, "w");

	for (const double x : vec)
	{
		fprintf(f, "%8.4lf\n", x);
	}
	
	fclose(f);
}
