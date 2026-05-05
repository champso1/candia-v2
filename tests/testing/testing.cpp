#include "Candia-v2/Common.hpp"
#include "Candia-v2/SplittingFn.hpp"
using namespace Candia2;

int main()
{
	getLogOptions().verbosity = LOG_INFO;
	
	auto compute_diff = [](double x, double y){ return 100.0*std::abs((x-y)/((x+y)/2.0));};
	SplittingFunction::setN3LOApproxType(3);
	SplittingFunction::update(4, 0, 0.0);
    
	std::unique_ptr<SplittingFunction> p3nsm_approx = std::make_unique<P3nsm>();
	std::unique_ptr<SplittingFunction> p3nsp_approx = std::make_unique<P3nsp>();
	std::unique_ptr<SplittingFunction> p3nsv_approx = std::make_unique<P3nsv>();
    
	std::unique_ptr<SplittingFunction> p3nsm_exact = std::make_unique<p3_exact::P3nsm>();
	std::unique_ptr<SplittingFunction> p3nsp_exact = std::make_unique<p3_exact::P3nsp>();
	std::unique_ptr<SplittingFunction> p3nsv_exact = std::make_unique<p3_exact::P3nsv>();
	
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
	for (double x = 1e-5; x<1.0; x+= (1.0-1e-5)/200.0) {
	 	double diff = compute_diff(p3nsv_approx->calcRegular(x), p3nsv_exact->calcRegular(x));
		std::cout << "x=" << x << "  ->  " << diff << "%\n";
	}
	// double diff = compute_diff(p3nsp_approx->calcDelta(), p3nsp_exact->calcDelta());
	// std::cout << "  ->  " << diff << "%\n";
}
	
