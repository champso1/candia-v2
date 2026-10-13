/**
 *  @file LHAPDFGrid.hpp
 *  @brief Contains the @a LHAPDFGrid class which runs the evolution many times to spit out an LHAPDF grid
 */

#pragma once

#include "Candia-v2/Candia.hpp"

#include <map>

namespace Candia2
{
	/**
	 *  @brief performs the evolution enough times to fill an lhapdf grid
	 *  the output is a pdf "set", consisting of one member only,
	 *  that is compatible with LHAPDF and able to be loaded via its API
	 *
	 *  Every evolution starts at the initial energy of the distribution, @a Distribution::Q0().
	 *  Each heavy-quark mass strictly inside the requested range of energies starts a new subgrid,
	 *  so the grid keeps the discontinuity of the PDFs (and of alpha_s) at the threshold:
	 *  the threshold appears as the last Q knot of the lower subgrid (pre-matching)
	 *  and as the first Q knot of the upper subgrid (post-matching).
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

		/** @brief one Q knot of the grid */
		struct Knot
		{
			double q{};      //!< energy
			uint nf{};       //!< number of active flavours used for this knot
			double alpha_s{}; //!< value of alpha_s at @a q (with @a nf flavours)
			std::map<int, ArrayGrid> pdfs{}; //!< x*f(x) on the grid, keyed by PDG id
		};
		std::vector<std::vector<Knot>> _subgrids{}; //!< one list of knots per flavour-number region, ordered in Q
		std::vector<double> _xvals{};

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
		/**
		 *  @param name name of the PDF set (and of the directory it is written to)
		 *  @param infofile_in_path optional template for the info file, providing metadata such as
		 *         SetDesc, Authors or Reference. Keys that LHAPDFGrid generates from the grid
		 *         (flavours, orders, ranges, masses, alpha_s) are ignored. Pass an empty path for none.
		 *  @param dist initial distribution; every evolution starts at @a dist.Q0()
		 */
		LHAPDFGrid(
			std::string const& name, std::filesystem::path const& infofile_in_path,
			Distribution& dist, Grid const& grid,
		    uint order, uint iterations, uint trunc_idx, double mur2_muf2)
			: _name{name}, _infofile_in_path{infofile_in_path},
			  _order{order}, _iterations{iterations}, _trunc_idx{trunc_idx}, _mur2_muf2{mur2_muf2},
			  _dist{dist}, _grid{grid}
		{
			if (!infofile_in_path.empty() && !std::filesystem::exists(infofile_in_path))
				log(LOG_ERROR, "LHAPDFGrid", "infofile_in_path is invalid ({})", infofile_in_path.string());
		}

		/** @brief evolves to every energy in @a qvals (which must not be below @a Distribution::Q0()) */
		inline void evolve(
		    std::vector<double> const& qvals,
			DGLAPSolver::options_type const& dglap_options)
		{
			_evolve_function(qvals, dglap_options, EvolType::Exact);
		}
		/** @brief as @a evolve(), with the truncated ansatz */
		inline void evolveTrunc(
		    std::vector<double> const& qvals,
			DGLAPSolver::options_type const& dglap_options)
		{
			_evolve_function(qvals, dglap_options, EvolType::Truncated);
		}
		/** @brief writes the info file and the grid of the last evolution to the directory @a name */
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

		/** @brief the PDG ids written to the grid, in order */
		std::vector<int> pids() const;
	};
}
