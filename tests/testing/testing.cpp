#include "arraygrid2.hpp"
#include "Candia-v2/ArrayGrid.hpp"
using namespace Candia2;

#include <random>
#include <charconv>
#include <cstring>

static constexpr uint ITERS = 15;
static constexpr uint GRID = 200;

void main1();
void main2() {}

int main(int argc, char** argv)
{
	if (argc != 2)
		return 1;
	int type{};
	std::from_chars(argv[1], argv[1]+strlen(argv[1]), type);
	switch (type) {
		case 1: main1(); break;
		case 2: main2(); break;
		default: throw std::runtime_error("unreachable");
	}
}

void main1()
{
	getLogOptions().show_debug_messages = true;
	std::random_device dev{};
	std::mt19937 mt{dev()};
	std::uniform_real_distribution<double> dist{};
	auto get_rand = [&](){ return dist(mt); };

	ArrayGrid2<6> arrgrid2({DISTS, 2, ITERS, ITERS, ITERS, GRID});

	for (uint j=0; j<DISTS; ++j) {
		for (uint s=0; s<2; ++s) {
			for (uint t=0; t<ITERS; ++t) {
				for (uint m=0; m<ITERS; ++m) {
					for (uint n=0; n<ITERS; ++n) {
						for (uint k=0; k<GRID; ++k)
							arrgrid2(j,s,t,m,n,k) = get_rand();
					}
				}
			}
		}
	}
}

