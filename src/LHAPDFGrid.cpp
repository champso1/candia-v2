// LHAPDFGrid.cpp

#include "Candia-v2/LHAPDFGrid.hpp"
#include "Candia-v2/Common.hpp"

#include <algorithm>

namespace Candia2
{
	void LHAPDFGrid::evolve(
		std::vector<double> const& qvals,
		DGLAPSolver::options_type const& dglap_options)
	{
		double q0 = qvals.front();
		double qf = qvals.back();
		AlphaS alphas_all(_order, q0, qf, _dist.alpha0(), _mur2_muf2);
		alphas_all.setVFNS(_dist.masses(), _dist.nfi(), _dist.nff());
		std::vector<std::pair<double,double>> as_qvals = alphas_all.getValues(qvals);
		std::vector<double> as_vals(as_qvals.size());
		std::ranges::transform(
			as_qvals, std::ranges::begin(as_vals),
			[](std::pair<double,double> const& p) -> double { return p.second; });
		_as_qs = std::move(qvals);
		_as_vals = std::move(as_vals);

		_xvals = _grid.points();

		log(LOG_DEBUG, "DGLAPSolverLHAPDF", "Energy values: {}", vec_to_str(_as_qs));

		auto enumerate =
			std::views::iota(uint{0}, _as_qs.size())
			| std::views::transform(
				[&](uint i){
					return std::make_pair(i, _as_qs[i]);
				});

		startLogIterations();
		for (auto[i, q] : enumerate) {
			getLogOptions().verbosity = LOG_INFO;
			logIterations(i, _as_qs.size()-1, "Evolutions");
			getLogOptions().verbosity = LOG_WARNING;
			
			_dist.setup(q0, q);
			AlphaS alphas(_order, q0, q, _dist.alpha0(), _mur2_muf2);
			alphas.setVFNS(_dist.masses(), _dist.nfi(), _dist.nff());
			// alphas.setFFNS(4);

			DGLAPSolver solver(_order, _grid, alphas, q, _iterations, _trunc_idx, _dist, _mur2_muf2);
			solver.getOptions() = dglap_options;
			std::vector<ArrayGrid> F = solver.evolve();

			std::unordered_map<int, ArrayGrid> map{};
			auto map_emplacer = [&](int k, ArrayGridView v) {
				auto emplace_res = map.emplace(k, ArrayGrid(v.size()));
				if (!emplace_res.second)
					log(LOG_ERROR, "LHAPDFGrid::evolve()", "map emplace failed");
				auto it = emplace_res.first;
				std::ranges::copy(v, it->second.begin());
			};

			map_emplacer(-5, F[11].view());
			map_emplacer(-4, F[10].view());
			map_emplacer(-3, F[9].view());
			map_emplacer(-2, F[8].view());
			map_emplacer(-1, F[7].view());
			map_emplacer(1, F[1].view());
			map_emplacer(2, F[2].view());
			map_emplacer(3, F[3].view());
			map_emplacer(4, F[4].view());
			map_emplacer(5, F[5].view());
			map_emplacer(21, F[0].view());

			std::vector<ArrayGrid> subtraction_pdfs = solver.calculateSubtractionPDFs();
			if (!subtraction_pdfs.empty()) {
				map_emplacer(FTILDE1, subtraction_pdfs[0].view());
				map_emplacer(FTILDE2, subtraction_pdfs[1].view());
				map_emplacer(FTILDE3, subtraction_pdfs[2].view());
				map_emplacer(FTILDENNLO, subtraction_pdfs[3].view());
				map_emplacer(FTILDEN3LO, subtraction_pdfs[4].view());
				map_emplacer(DELTAF1, subtraction_pdfs[5].view());
				map_emplacer(DELTAF2, subtraction_pdfs[6].view());
				map_emplacer(DELTAF3, subtraction_pdfs[7].view());
			}

		    _all_pdfs.emplace_back(std::make_pair(q, map));
		}

		getLogOptions().verbosity = LOG_INFO;
		endLogIterations();
		getLogOptions().verbosity = LOG_WARNING;
	}
	
	void LHAPDFGrid::evolve(
		double q0, double qf, double dq,
		DGLAPSolver::options_type const& dglap_options)
	{
	    double nsize = (qf-q0)/dq;
		uint size = std::trunc(nsize) + 1;
		std::vector<double> qvals(size);
		auto as_qs_view =
			std::views::iota(uint{0}, size)
			| std::views::transform([q0,dq](uint x) -> double { return q0 + dq*x; });
	    std::vector<double> as_qs(as_qs_view.begin(), as_qs_view.end());
	    if (std::ranges::find(as_qs, qf) == as_qs.end())
			as_qs.emplace_back(qf);
		evolve(as_qs, dglap_options);
	}

	void LHAPDFGrid::write()
	{
		log(LOG_DEBUG, "DGLAPSolverLHAPDF", "Spitting out the stuff...");
		std::filesystem::path pdfdir_path = std::filesystem::current_path()/_name;
		if (!std::filesystem::exists(pdfdir_path)) {
			if (!std::filesystem::create_directory(pdfdir_path))
				log(LOG_ERROR, "DGLAPSolverLHAPDF", "Failed to create the pdf outpout directory: '{}'", pdfdir_path.string());
		}
		
		std::string infofile_name = std::format("{}.info", _name);
		std::ifstream infofile_in(_infofile_in_path);
		std::string infofile_in_contents(
			std::istreambuf_iterator<char>(infofile_in),
			std::istreambuf_iterator<char>{});
		std::string pids_replace_str = "%PIDS%";
		std::string order_replace_str = "%ORDER%";
		std::string as_qs_replace_str = "%AS_QS%";
		std::string as_vals_replace_str = "%AS_VALS%";

		auto perform_replace = [](
			std::string& str,
			std::string const& what_to_replace,
			std::string const& replace_with)
		{
			std::string::size_type ipos;
			while ((ipos = str.find(what_to_replace)) != std::string::npos)
				str.replace(ipos, what_to_replace.size(), replace_with);
		};

		auto pids =
			_all_pdfs[0].second
			| std::views::transform(
				[](std::pair<const int, ArrayGrid> const& p) -> int {
					return p.first;
				});

		perform_replace(infofile_in_contents, pids_replace_str, vec_to_str2(pids));
		perform_replace(infofile_in_contents, order_replace_str, std::to_string(_order));
		perform_replace(infofile_in_contents, as_qs_replace_str, vec_to_str2(_as_qs));
		perform_replace(infofile_in_contents, as_vals_replace_str, vec_to_str2(_as_vals));

		std::filesystem::path infofile_path = pdfdir_path/infofile_name;
		std::ofstream infofile(infofile_path);
		std::copy(infofile_in_contents.begin(), infofile_in_contents.end(), std::ostreambuf_iterator<char>(infofile));
		
		std::string datafile_name = std::format("{}_0000.dat", _name);
		std::filesystem::path datafile_path = pdfdir_path/datafile_name;
		std::ofstream datafile(datafile_path);
		datafile << "PdfType: central\n";
		datafile << "Format: lhagrid1\n";
		datafile << "---\n";
		datafile << "   " << vec_to_str2(_xvals, " ") << '\n';
		datafile << "   " << vec_to_str2(_as_qs, " ") << '\n';
		datafile << " " << vec_to_str2(pids, " ") << '\n';

		// TODO: what the fuck is this nonsense
	    // std::vector<std::pair<double, std::map<int,ArrayGrid>>>
		datafile << std::setprecision(10) << std::scientific;
		for (uint ix=0; ix<_xvals.size(); ++ix) {
			for (uint iq=0; iq<_as_qs.size(); ++iq) {
				datafile << "   ";
				for (int pid : pids)
					datafile << std::setw(17) << _all_pdfs[iq].second[pid](ix) << ' ';
				datafile << '\n';
			}
		}
		datafile << "---\n";
	}
}
