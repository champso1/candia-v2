#include "Candia-v2/Couplings.hpp"

namespace Candia2
{
	std::array<AlphaQED::nfud_threshold_type,9> AlphaQED::_nfud_thresholds
	{
		nfud_threshold_type{0,0,0},
		nfud_threshold_type{1,0,0},
		nfud_threshold_type{2,0,0},
		nfud_threshold_type{2,2,1},
		nfud_threshold_type{2,2,2},
		nfud_threshold_type{3,2,2},
		nfud_threshold_type{3,3,2},
		nfud_threshold_type{3,3,3},
	};
}
