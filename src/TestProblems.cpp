#include "TestProblems.hpp"

using namespace TestProblems;

template <size_t dim>
Problem<dim> create_problem(ProblemName p)
{
	if (p == TOWNSEND && dim == 2)
	{
		auto f = [](const RealVector<dim> &x)
		{
			double ret = 0;
			ret -= (cos((x[0] - 0.1) * x[1]) * cos((x[0] - 0.1) * x[1]));
			ret -= x[0] * sin(3 * x[0] + x[1]);
			return ret;
		};

		RealVector<dim> lb(-2.25, -2.5);
		RealVector<dim> ub(2.25, 1.75);

		Problem<dim> problem(f, lb, ub);

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double t = atan2(x[0], x[1]);
			double ret = 0;
			ret += x[0] * x[0] + x[1] * x[1];
			ret -= (2 * cos(t) - cos(2*t)/2 - cos(3*t)/4 - cos(4*t)/8) * (2 * cos(t) - cos(2*t)/2 - cos(3*t)/4 - cos(4*t)/8);
			ret -= (2 * sin(t)) * (2 * sin(t));
			return ret; });

		return problem;
	}
	else
		throw std::runtime_error("Problem not implemented.\n \
			Check the documentation for the available problems.\n \
			Check if the dim template parameter is consistent with the specified problem.");
}

template <size_t dim>
double get_exact_value(ProblemName p)
{
	if (p == TOWNSEND && dim == 2)
	{
		return -2.0239884;
	}
	else
		throw std::runtime_error("Problem not implemented.\n \
			Check the documentation for the available problems.\n \
			Check if the dim template parameter is consistent with the specified problem.");
}

template <size_t dim>
RealVector<dim> get_exact_position(ProblemName p)
{
	if (p == TOWNSEND && dim == 2)
	{
		return RealVector<dim>(-3.1302468, -1.5821422);
	}
	else
		throw std::runtime_error("Problem not implemented.\n \
			Check the documentation for the available problems.\n \
			Check if the dim template parameter is consistent with the specified problem.");
}
