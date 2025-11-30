#include <print>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <fstream>
#include <ostream>
#include <iomanip>
using namespace std;
using uint = unsigned;

const static vector<double> xtab{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
const static uint nx = 1000;

pair<vector<double>,vector<int>> gen_grid0()
{
	const uint xtab_len = xtab.size();
	vector<double> Ntab(xtab_len);
	vector<int> ntab(xtab_len);

	double temp = -log10(xtab[0]);

	for (uint i=1; i<xtab_len; i++)
		Ntab[i] = (double)(nx-1)*log10(xtab[i]/xtab[i-1])/temp;

	ntab[0] = nx-1;

	for (uint i=1; i<xtab_len; i++)
	{
		ntab[i] =  (int)Ntab[i];
		Ntab[i] -= (double)ntab[i];
		ntab[0] -= ntab[i];
	}

	for (uint i=1; i<xtab_len; i++)
	{
		if (ntab[i] == 0)
		{
			ntab[i] =  1;
			Ntab[i] -= 1.0;
			ntab[0] -= 1;
		}
	}

	uint n;
	for ( ; ntab[0]<0; ntab[0]++)
	{
		n=0;

		for (uint i=1; i<xtab_len; i++)
		{
			if (ntab[i] != 1)
			{
				if ((n == 0) || (Ntab[i] <= temp))
				{
					n = i;
					temp = Ntab[i];
				}
			}
		}

		ntab[n]--;
		Ntab[n] += 1.0;
	}

	for ( ; ntab[0]>0; ntab[0]--)
	{
		n=0;

		for (uint i=1; i<xtab_len; i++)
		{
			if ((n == 0) || (Ntab[i] > temp))
			{
				n = i;
				temp = Ntab[i];
			}
		}

		ntab[n]++;
		Ntab[n] -= 1.0;
	}

	for (uint i=1; i<xtab_len; i++)
		ntab[i] += ntab[i-1];

	vector<double> points(nx);
	double lstep;
	for (uint i=0; i<xtab_len-1; i++)
	{
		lstep=log10(xtab[i+1]/xtab[i])/(double)(ntab[i+1]-ntab[i]);

		for (int j=ntab[i]; j<ntab[i+1]; j++)
			points.at(j) = xtab[i]*pow(10.0, lstep*(double)(j-ntab[i]));
	}

	points.at(nx-1) = 1.0;

	return {points, ntab};
}

pair<vector<double>,vector<int>> gen_grid3()
{
	vector<double> x(nx-xtab.size()+1);
	const double xmin = 1.0e-7;
	const double xmax = 1.0;
	const double ymin = log10(xmin);
	const double ymax = log10(xmax);
	auto func = [=](double u) -> double
	{
		double fac1 = 0.25*(0.5-u)*(0.5-u);
		return ymin + (ymax-ymin)*fac1;
	};
	double func_min = func(0.5);
	double func_max = func(1.0);
	for (uint i=0; i<nx-xtab.size()+1; ++i)
	{
		const double u = static_cast<double>(i)/(static_cast<double>(nx) - 1.0);
		// double _y = ymin + (ymax-ymin)*0.5*(1.0 - cos(PI*u));
		double _y = ymin + (ymax-ymin)/(func_max - func_min)*func(u);
		x[i] = pow(10.0, _y);
	}

	// insert the tabulated points
	ranges::sort(x);
	set<double> set(x.begin(), x.end());
	set.insert(xtab.begin(), xtab.end());
	x = vector<double>(set.begin(), set.end());

	// build the ntab array
	vector<int> ntab{};
	for (const double _x : xtab)
	{
		auto it = ranges::lower_bound(x, _x);
		if (it != x.end() && abs(*it - _x) < 1e-14)
			ntab.emplace_back(distance(x.begin(), it));
	}

	// fix front and back 
	x.front() = xmin;
	if (x.back() != 1.0)
		x.back() = 1.0;

	return {x, ntab};
}

pair<vector<double>,vector<int>> gen_grid4()
{
    vector<double> x(nx-xtab.size()+1);

	uint mid = x.size()*.75;
	println("Total size of array is given by {}", x.size());
	println("\"Midpoint\" for the values is {}", mid);
	
	double xmin = xtab.front();
	double xmax = 0.5;
	double ymin = log10(xmin);
	double ymax = log10(xmax);
	println("[{},{}] -> [{},{}]", xmin, xmax, ymin, ymax);
	
	for (uint i=0; i<mid; ++i)
	{
		double u = static_cast<double>(i)/(static_cast<double>(mid)-1.0);
		double _y = ymin + (ymax-ymin)*u;
		x[i] = pow(10.0, _y);

		println("{}: u={:.5} y={:.5} x={:.5}", i, u, _y, x[i]);
	}

	println("-------------------- SWITCHING --------------------");

    xmin = 0.5;
    xmax = xtab.back();
    ymin = log10(xmin);
    ymax = log10(xmax);
	println("[{},{}] -> [{},{}]", xmin, xmax, ymin, ymax);
	
	for (uint i=mid; i<x.size(); ++i)
	{
		double u = static_cast<double>(i)/(static_cast<double>(x.size()-mid)-1.0);
		double _y = ymin + (ymax-ymin)*u;
		x[i] = pow(10.0, _y);

		println("{}: u={:.5} y={:.5} x={:.5}", i, u, _y, x[i]);
	}

	// insert the tabulated points
	ranges::sort(x);
	set<double> set(x.begin(), x.end());
	set.insert(xtab.begin(), xtab.end());
	x = vector<double>(set.begin(), set.end());

	// build the ntab array
	vector<int> ntab{};
	for (const double _x : xtab)
	{
		auto it = ranges::lower_bound(x, _x);
		if (it != x.end() && abs(*it - _x) < 1e-14)
			ntab.emplace_back(distance(x.begin(), it));
	}

	// fix front and back 
	x.front() = xmin;
	if (x.back() != 1.0)
		x.back() = 1.0;

	return {x, ntab};
}

void print_vals(vector<int> ntab, vector<double> x)
{
	println("Ntab values:");
    ranges::for_each(ntab, [](int x){ print("{} ", x); });
	println();
	println("Xtab values:");
    ranges::for_each(ntab, [&X=x](int _x){ print("{} ", X[_x]); });
	println();
	println("X values:");
	ranges::for_each(x, [](double x){ println("  {}", x); });
	println();
}

void save_vals(vector<double> x)
{
	ofstream outfile{"points.dat"};
	outfile << fixed << setprecision(10);
	ranges::for_each(x, [&os = outfile](const double x) { println(os, "{:.10e}", x); });
}


int main()
{
	// auto [x1, ntab1] = gen_grid0();
	// auto [x3, ntab3] = gen_grid3();
	auto [x4, ntab4] = gen_grid4();

	// print_vals(ntab1, x1);
	// print_vals(ntab3, x3);
	print_vals(ntab4, x4);

	// save_vals(x1);
	// save_vals(x3);
	// save_vals(x4);
}

