/**
 *  @file LHAPDFGrid.hpp
 *  @brief Contains the @a LHAPDFGrid class which runs the evolution many times to spit out an LHAPDF grid
 */

#pragma once

#include "Candia-v2/Candia.hpp"

namespace Candia2
{
	/**
	 *  @brief performs the evolution enough times to fill an lhapdf grid
	 *  the output is a pdf "set", consisting of one file only,
	 *  that is compatible with LHAPDF and able to be loaded via its API
	 */
	class LHAPDFGrid final
	{
		std::string _name;
		std::filesystem::path _infofile_in_path;

	    uint _order{3};
	    uint _iterations{10};
	    uint _trunc_idx{10};
	    double _mur2_muf2{1.0};

		Distribution& _dist;
		Grid _grid;

	    std::vector<std::pair<double, std::unordered_map<int,ArrayGrid>>> _all_pdfs{};
		std::vector<double> _as_qs, _as_vals, _xvals;

	public:
		enum SubtractionPDFIndices : int
		{
		    C_FTILDE1=0,
			C_FTILDE2   ,
			C_FTILDE3   ,
			C_FTILDENNLO, // ft1 + ft2
			C_FTILDEN3LO, // ft1 + ft2 + ft3
			C_DELTAF1   , // fc - ft1
			C_DELTAF2   , // fc - ft2
			C_DELTAF3   , // fc - ft3
			C_DELTAFNNLO, // fc - ftNNLO
			C_DELTAFN3LO, // fc - ftN3LO
			B_FTILDE1   ,
			B_FTILDE2   ,
			B_FTILDE3   ,
			B_FTILDENNLO, // ft1 + ft2
			B_FTILDEN3LO, // ft1 + ft2 + ft3
			B_DELTAF1   , // fb - ft1
			B_DELTAF2   , // fb - ft2
			B_DELTAF3   , // fb - ft3
			B_DELTAFNNLO, // fb - ftNNLO
			B_DELTAFN3LO  // fb - ftN3LO
		};
	private:
		uint _lhapdf_subpdf_pid_offset{9000};
	public:
		LHAPDFGrid(
			std::string const& name, std::filesystem::path const& infofile_in_path,
			Distribution& dist, Grid const& grid,
		    uint order, uint iterations, uint trunc_idx, double mur2_muf2)
			: _name{name}, _infofile_in_path{infofile_in_path},
			  _order{order}, _iterations{iterations}, _trunc_idx{trunc_idx}, _mur2_muf2{mur2_muf2},
			  _dist{dist}, _grid{grid}
		{
			if (!std::filesystem::exists(infofile_in_path))
				log(LOG_ERROR, "DGLAPSolverLHAPDF", "infofile_in_path is invalid ({})", infofile_in_path.string());
		}
		
		inline void evolve(
		    std::vector<double> const& qvals,
			DGLAPSolver::options_type const& dglap_options)
		{
			_evolve_function(qvals, dglap_options, EvolType::Exact);
		}
		inline void evolveTrunc(
		    std::vector<double> const& qvals,
			DGLAPSolver::options_type const& dglap_options)
		{
			_evolve_function(qvals, dglap_options, EvolType::Truncated);
		}
		void write();
	private:
		enum class EvolType
		{
			Exact,
			Truncated,
		};
		
		void _evolve_function(
			std::vector<double> const& qvals,
			DGLAPSolver::options_type const& dglap_options,
			EvolType evol_type);
	};
}
