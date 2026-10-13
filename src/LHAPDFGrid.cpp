// LHAPDFGrid.cpp

#include "Candia-v2/LHAPDFGrid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>
#include <format>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

namespace Candia2
{
	namespace
	{
		/** @brief restores the global log verbosity when it goes out of scope */
		struct VerbosityGuard final
		{
			int saved{getLogOptions().verbosity};
			~VerbosityGuard() { getLogOptions().verbosity = saved; }
		};

		/** @brief a Q knot to evolve to, and the number of flavours to use there */
		struct KnotSpec final
		{
			double q{};
			uint nf{};
		};
	}

	void LHAPDFGrid::_evolve_function(
		std::vector<double> const& qvals,
		DGLAPSolver::options_type const& dglap_options,
		LHAPDFGrid::EvolType evol_type)
	{
		if (qvals.empty())
			log(LOG_ERROR, "LHAPDFGrid", "No energies to evolve to");

		std::vector<double> knots{qvals};
		std::ranges::sort(knots);
		auto duplicates = std::ranges::unique(knots);
		knots.erase(duplicates.begin(), duplicates.end());

		// every evolution starts at the initial energy of the distribution
		const double q0 = _dist.Q0();
		const double dist_qf = _dist.Qf();
		const double qmin = knots.front();
		const double qmax = knots.back();
		if (qmin < q0)
			log(LOG_ERROR, "LHAPDFGrid", "Cannot evolve to Q={} below the initial energy of the distribution (Q0={})", qmin, q0);

		_dist.setup(q0, qmax);
		const uint nfi = _dist.nfi();
		auto const masses = _dist.masses();

		// heavy-quark thresholds above the initial flavours:
		// those at or below the first knot have already been crossed there,
		// those strictly inside the range of knots start a new subgrid
		uint nf_start = nfi;
		std::vector<double> boundaries{};
		for (uint h=nfi+1; h<=6; ++h) {
			const double m = masses[h];
			if (m <= 0.0 || m >= qmax)
				break;
			if (m <= qmin)
				++nf_start;
			else
				boundaries.push_back(m);
		}
		const uint nf_max = nf_start + static_cast<uint>(boundaries.size());

		// subgrid i has nf_start+i flavours and runs from boundaries[i-1] to boundaries[i];
		// each threshold is the last knot of one subgrid and the first knot of the next
		std::vector<std::vector<KnotSpec>> specs(boundaries.size()+1);
		constexpr double inf = std::numeric_limits<double>::infinity();
		uint num_evolutions = 0;
		for (uint i=0; i<specs.size(); ++i) {
			const uint nf = nf_start + i;
			const double lo = (i == 0) ? -inf : boundaries[i-1];
			const double hi = (i == boundaries.size()) ? inf : boundaries[i];
			if (i > 0)
				specs[i].push_back({lo, nf});
			for (double q : knots) {
				if (q > lo && q < hi)
					specs[i].push_back({q, nf});
			}
			if (i < boundaries.size())
				specs[i].push_back({hi, nf});

			if (specs[i].size() < 4) {
				log(LOG_WARNING, "LHAPDFGrid",
					"The subgrid with {} flavours has only {} Q knot(s); LHAPDF interpolates more accurately with at least 4",
					nf, specs[i].size());
			}
			num_evolutions += static_cast<uint>(specs[i].size());
		}

		log(LOG_DEBUG, "LHAPDFGrid", "Q0 = {}, flavours {} to {}, thresholds inside the grid: {}",
			q0, nf_start, nf_max, vec_to_str2(boundaries));

		_subgrids.clear();
		_xvals = _grid.points();

		// the solvers print a lot at LOG_INFO, so only the progress bar is shown at that level
		VerbosityGuard verbosity_guard{};
		const int quiet = (verbosity_guard.saved == LOG_NONE) ? LOG_NONE : std::max(verbosity_guard.saved, static_cast<int>(LOG_WARNING));

		uint done = 0;
		startLogIterations();
		for (auto const& subgrid_specs : specs) {
			auto& subgrid = _subgrids.emplace_back();
			for (auto [q, nf] : subgrid_specs) {
				getLogOptions().verbosity = quiet;

				_dist.setup(q0, q);
				AlphaS alphas(_order, q0, q, _dist.alpha0(), _mur2_muf2);
				// the number of flavours is set explicitly so that a threshold knot can be
				// evaluated both just below (pre-matching) and just above (post-matching)
				alphas.setVFNS(_dist.masses(), _dist.nfi(), nf);

				DGLAPSolver solver(_order, _grid, alphas, q, _iterations, _trunc_idx, _dist, _mur2_muf2);
				solver.getOptions() = dglap_options;
				auto const& F = (evol_type == EvolType::Exact) ? solver.evolve() : solver.evolveTrunc();

				// alpha_s at q: with nff = nf, masses[nf+1] is set to q, so this is its "pre-threshold" value
				Knot& knot = subgrid.emplace_back(Knot{q, nf, alphas.pre(nf+1), {}});
				auto add = [&knot](int pid, ArrayGrid const& values) { knot.pdfs.emplace(pid, values); };

				// PDG ids: 1 = d, 2 = u, ...; internally index DIST_U = 1 and DIST_D = 2, antiquarks at +6
				add(21, F[0]);
				for (auto [pid, idx] : {std::pair{1, DIST_D}, {2, DIST_U}, {3, DIST_S}, {4, DIST_C}, {5, DIST_B}}) {
					add(pid, F[idx]);
					add(-pid, F[idx+6]);
				}
				if (nf_max == 6) {
					add(6, F[DIST_T]);
					add(-6, F[DIST_T+6]);
				}

				std::vector<ArrayGrid> subtraction_pdfs = solver.calculateSubtractionPDFs();
				for (uint i=0; i<subtraction_pdfs.size(); ++i)
					add(static_cast<int>(_lhapdf_subpdf_pid_offset + i), subtraction_pdfs[i]);

				getLogOptions().verbosity = verbosity_guard.saved;
				logIterations(++done, num_evolutions, "LHAPDFGrid");
			}
		}
		getLogOptions().verbosity = verbosity_guard.saved;
		endLogIterations();

		// leave the distribution as it was given to us
		_dist.setup(q0, dist_qf);
	}

	std::vector<int> LHAPDFGrid::pids() const
	{
		std::vector<int> out{};
		if (_subgrids.empty() || _subgrids.front().empty())
			return out;
		for (auto const& [pid, values] : _subgrids.front().front().pdfs)
			out.push_back(pid);
		return out;
	}

	void LHAPDFGrid::write()
	{
		if (_subgrids.empty())
			log(LOG_ERROR, "LHAPDFGrid", "Nothing to write: call evolve() or evolveTrunc() first");

		namespace fs = std::filesystem;
		fs::path pdfdir_path = fs::current_path()/_name;
		std::error_code ec{};
		fs::create_directories(pdfdir_path, ec);
		if (ec)
			log(LOG_ERROR, "LHAPDFGrid", "Failed to create the PDF output directory '{}': {}", pdfdir_path.string(), ec.message());

		const std::vector<int> pid_list = pids();
		auto const& masses = _dist.masses();

		std::vector<double> as_qs{}, as_vals{};
		uint nf_max = 0;
		for (auto const& subgrid : _subgrids) {
			for (auto const& knot : subgrid) {
				as_qs.push_back(knot.q);
				as_vals.push_back(knot.alpha_s);
				nf_max = std::max(nf_max, knot.nf);
			}
		}

		// ---------- info file ----------
		auto list = [](auto const& values) { return std::format("[{}]", vec_to_str2(values)); };
		auto number = [](double v) { return std::format("{:.10e}", v); };
		const std::vector<std::pair<std::string, std::string>> generated{
			{"Format", "lhagrid1"},
			{"Flavors", list(pid_list)},
			{"OrderQCD", std::to_string(_order)},
			{"FlavorScheme", "variable"},
			{"NumFlavors", std::to_string(nf_max)},
			{"XMin", number(_xvals.front())},
			{"XMax", number(_xvals.back())},
			{"QMin", number(as_qs.front())},
			{"QMax", number(as_qs.back())},
			{"MUp", number(masses[DIST_U])},
			{"MDown", number(masses[DIST_D])},
			{"MStrange", number(masses[DIST_S])},
			{"MCharm", number(masses[DIST_C])},
			{"MBottom", number(masses[DIST_B])},
			{"MTop", number(masses[DIST_T])},
			{"AlphaS_OrderQCD", std::to_string(_order)},
			{"AlphaS_Type", "ipol"},
			{"AlphaS_Qs", list(as_qs)},
			{"AlphaS_Vals", list(as_vals)},
		};
		const std::vector<std::pair<std::string, std::string>> defaults{
			{"SetDesc", std::format("'Candia-v2 evolution ({})'", _name)},
			{"DataVersion", "1"},
			{"NumMembers", "1"},
		};
		std::set<std::string> generated_keys{};
		for (auto const& [key, value] : generated)
			generated_keys.insert(key);

		// the template provides metadata; any key we generate is dropped from it,
		// together with its continuation lines (indented, or list items)
		std::string info_text{};
		std::set<std::string> template_keys{};
		if (!_infofile_in_path.empty()) {
			std::istringstream template_in(read_file(_infofile_in_path));
			std::string line{};
			bool dropping = false;
			while (std::getline(template_in, line)) {
				if (!line.empty() && line.back() == '\r')
					line.pop_back();
				const bool continuation = !line.empty() && (line[0] == ' ' || line[0] == '\t' || line[0] == '-');
				if (!continuation) {
					const auto colon = line.find(':');
					const std::string key = (colon == std::string::npos) ? std::string{} : line.substr(0, colon);
					dropping = generated_keys.contains(key);
					if (!key.empty() && !dropping)
						template_keys.insert(key);
				}
				if (!dropping)
					info_text += line + '\n';
			}
		}
		for (auto const& [key, value] : defaults) {
			if (!template_keys.contains(key))
				info_text += std::format("{}: {}\n", key, value);
		}
		for (auto const& [key, value] : generated)
			info_text += std::format("{}: {}\n", key, value);

		fs::path infofile_path = pdfdir_path/std::format("{}.info", _name);
		std::ofstream infofile(infofile_path);
		infofile << info_text;
		if (!infofile)
			log(LOG_ERROR, "LHAPDFGrid", "Failed to write '{}'", infofile_path.string());

		// ---------- grid file: one block per subgrid ----------
		fs::path datafile_path = pdfdir_path/std::format("{}_0000.dat", _name);
		std::ofstream datafile(datafile_path);
		datafile << "PdfType: central\n";
		datafile << "Format: lhagrid1\n";
		datafile << "---\n";
		datafile << std::setprecision(10) << std::scientific;
		for (auto const& subgrid : _subgrids) {
			std::vector<double> qs{};
			for (auto const& knot : subgrid)
				qs.push_back(knot.q);

			datafile << "   " << vec_to_str2(_xvals, " ") << '\n';
			datafile << "   " << vec_to_str2(qs, " ") << '\n';
			datafile << " " << vec_to_str2(pid_list, " ") << '\n';
			// x is the outer index, Q the inner one
			for (uint ix=0; ix<_xvals.size(); ++ix) {
				for (auto const& knot : subgrid) {
					datafile << "   ";
					for (int pid : pid_list)
						datafile << std::setw(17) << knot.pdfs.at(pid)[ix] << ' ';
					datafile << '\n';
				}
			}
			datafile << "---\n";
		}
		if (!datafile)
			log(LOG_ERROR, "LHAPDFGrid", "Failed to write '{}'", datafile_path.string());
	}
}
