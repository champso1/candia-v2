#ifndef __EXPRESSION_HPP
#define __EXPRESSION_HPP

namespace Candia2
{

	/** Class that handles expressions that contain regular parts,
	 *  singular parts (associated with plus distributions),
	 *  and "local" parts, associated with delta functions
	 */
	class Expression
	{
	protected:
		Expression() = default;
	public:
		virtual ~Expression() = default;
		
		/** @brief Regular part of the expression
		 */
		virtual double Regular(double x) const;

		/** @brief Singular (plus-distribution) part of the expression
		 */
		virtual double Plus(double x) const;

		/** @brief Delta function part of the distribution
		 */
		virtual double Delta(double x) const;
		
	};
	
}; // namespace Candia2

#endif // __EXPRESSION_HPP
