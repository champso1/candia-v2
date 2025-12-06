#ifndef __EXPRESSION_HPP
#define __EXPRESSION_HPP

namespace Candia2
{
	class Expression
	{
	protected:
		Expression() = default;
	public:
		virtual ~Expression() = default;

		virtual double regular(double x) const;
		virtual double plus(double x) const;
		virtual double delta(double x) const;
		
	};
	
};

#endif
