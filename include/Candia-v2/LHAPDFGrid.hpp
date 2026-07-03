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
		enum AdditionalSubtractionPDFs : int
		{
		    C_FTILDE1    = 9001,
			C_FTILDE2    = 9002,
			C_FTILDE3    = 9003,
			C_FTILDENNLO = 9004, // ft1 + ft2
			C_FTILDEN3LO = 9005, // ft1 + ft2 + ft3
			C_DELTAF1    = 9006, // fc - ft1
			C_DELTAF2    = 9007, // fc - ft2
			C_DELTAF3    = 9008, // fc - ft3
			C_DELTAFNNLO = 9009, // fc - ftNNLO
			C_DELTAFN3LO = 9010, // fc - ftN3LO
			B_FTILDE1    = 9101,
			B_FTILDE2    = 9102,
			B_FTILDE3    = 9103,
			B_FTILDENNLO = 9104, // ft1 + ft2
			B_FTILDEN3LO = 9105, // ft1 + ft2 + ft3
			B_DELTAF1    = 9106, // fb - ft1
			B_DELTAF2    = 9107, // fb - ft2
			B_DELTAF3    = 9108, // fb - ft3
			B_DELTAFNNLO = 9109, // fb - ftNNLO
			B_DELTAFN3LO = 9110  // fb - ftN3LO
		};
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
		
		void evolve(
			double q0, double qf, double dq,
			DGLAPSolver::options_type const& dglap_options);
		void evolve(
		    std::vector<double> const& qvals,
			DGLAPSolver::options_type const& dglap_options);
		void write();
	};
}
