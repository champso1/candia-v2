#include "Candia-v2/Common.hpp"
#include "Candia-v2/Couplings.hpp"
#include <iterator>
using namespace Candia2;

#include <vector>
#include <fstream>

int main()
{
	getLogOptions().verbosity = LOG_INFO;

	uint order = 0;
	double Q0 = 1e-2;
	double Qf = 1e2;
	double logq0 = std::log10(Q0);
	double logqf = std::log10(Qf);
	double mur2_muf2 = 1;
	double a0 = ALPHAQED_MTAU;
	uint nf = 4;
	uint nl = 3;
	AlphaQED qed(order, Q0, Qf, a0, mur2_muf2);
	qed.update(nf, nl);

	uint num_energies = 100;
	double num_energies_f = num_energies;
	double dq = (Qf-Q0)/(num_energies_f-1);
	double log_dq = (logqf-logq0)/(num_energies_f-1);
	
	auto values =
		std::views::iota(uint{0}, num_energies)
		| std::views::transform(
			[&](uint i){
				double q = std::pow(10.0, logq0 + static_cast<double>(i)*log_dq);
				double a = qed.evaluate(Q0, q, a0);
				return std::make_pair(q, a);
			});
	std::ofstream outfile("alphaqed.dat");
	for (auto [e,a] : values)
		outfile << e << ' ' << a << '\n';
}
