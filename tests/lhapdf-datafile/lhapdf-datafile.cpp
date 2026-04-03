#include "Candia-v2/Candia.hpp"
using namespace Candia2;
using namespace std;

int main()
{
	auto& log_options = getLogOptions();
	log_options.show_debug_messages = true;
	log_options.show_thread_output = true;
	
	std::unique_ptr<Distribution> dist = std::make_unique<LesHouchesDistribution>();
	std::filesystem::path infofile_in_path("infofile.in");
	DGLAPSolverLHAPDF solver(
		"testpdf", infofile_in_path,
		std::move(dist),
		0, 10, 10, 1.0);
	solver.evolve(std::numbers::sqrt2, 100, 10.0);
	solver.write();
}
