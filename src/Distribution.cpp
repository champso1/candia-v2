#include "Candia-v2/Distribution.hpp"

#include <cmath>

namespace Candia2
{
	
	
	double LesHouchesDistribution::xuv(const double x) const {
		return 5.1072*std::pow(x, 0.8)*std::pow(1.0-x, 3.0);
	}

	double LesHouchesDistribution::xdv(const double x) const {
		return 3.06432*std::pow(x, 0.8)*std::pow(1.0-x, 4.0);
	}

	double LesHouchesDistribution::xg(const double x) const {
		return 1.7*std::pow(x, -0.1)*std::pow(1.0-x, 5.0);
	}

	double LesHouchesDistribution::xdb(const double x) const {
		return 0.1939875*std::pow(x, -0.1)*std::pow(1.0-x, 6.0);
	}

	double LesHouchesDistribution::xub(const double x) const {
		return xdb(x)*(1.0-x);
	}

	double LesHouchesDistribution::xs(const double x) const {
		return 0.2*(xub(x) + xdb(x));
	}
	
}
