/** @file
 *
 *  Contains information related to initial PDF distributions.
 *  Only the Les Houche "toy" model has been implemented concretely,
 *  as it is what is used in current benchmarkings.
 */

#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP


namespace Candia2
{

	/** Represents an initial pdf distribution
	 */
	class Distribution
	{
	public:
		/** @name Constructors/destructors
		 */
		///@{
		Distribution() = default; //!< default constructor
		~Distribution() = default; //!< default destructor
		///@}

		/** @name Initial conditions
		 */
		///@{
		/** @brief up-valence */
		virtual double xuv(const double x) const = 0;

		/** @brief down-valence */
		virtual double xdv(const double x) const = 0;

		/** @brief gluon */
		virtual double xg (const double x) const = 0;

		/** @brief down-bar */
		virtual double xdb(const double x) const = 0;

	    /** @brief up-bar */
		virtual double xub(const double x) const = 0;

		/** @brief strange */
		virtual double xs (const double x) const = 0;
		///@}
	};



	/** LesHouches toy model initial distributions
	 */
	class LesHouchesDistribution : public Distribution
	{
		double xuv(const double x) const;
		double xdv(const double x) const;
		double xg (const double x) const;
		double xdb(const double x) const;
		double xub(const double x) const;
		double xs (const double x) const;
	};
	
}


#endif // __DISTRIBUTION_HPP
