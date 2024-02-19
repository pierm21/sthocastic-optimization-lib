#pragma once

#include "TestProblems.hpp"

template <size_t dim>
Problem<dim> TestProblems::create_problem(ProblemName p)
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
	else if (p == GOMEZ_LEVY && dim == 2)
	{
		auto f = [](const RealVector<dim> &x)
		{
			double ret = 0;
			double x2 = x[0] * x[0];
			double y2 = x[1] * x[1];
			ret += 4.0 * x2 - 2.1 * x2 * x2 + (1.0 / 3) * x2 * x2 * x2;
			ret += x[0] * x[1];
			ret += 4 * y2 * y2 - 4 * y2;
			return ret;
		};

		RealVector<dim> lb(-1, -1);
		RealVector<dim> ub(0.75, 1);

		Problem<dim> problem(f, lb, ub);

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = 0;
			ret -= std::sin(4 * M_PI * x[0]);
			ret += 2.0 * std::sin(2 * M_PI * x[1]) * std::sin(2 * M_PI * x[1]);
			ret -= 1.5;
			return ret; });

		return problem;
	}
	/*else if (p == _G10 && dim == 8)
	{
		auto f = [](const RealVector<dim> &x)
		{
			return x[0] + x[1] + x[2];
		};

		RealVector<dim> lb(100, 1000, 1000, 10, 10, 10, 10, 10);
		RealVector<dim> ub(10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000);

		Problem<dim> problem(f, lb, ub);

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = -1.0;
			ret += 0.0025 * (x[3] + x[5]);
			return ret; });

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = -1.0;
			ret += 0.0025 * (x[4] + x[6] - x[3]);
			return ret; });

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = -1.0;
			ret += 0.01 * (x[7] - x[4]);
			return ret; });

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = - x[0] * x[5];
			ret += 833.33252 * x[3] + 100 * x[0] - 83333.333;
			return ret; });

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = - x[1] * x[6];
			ret += 1250 * x[4] + x[1] * x[3] - 1250 * x[4];
			return ret; });

		problem.add_inequality_constraint([](const RealVector<dim> &x)
										  {
			double ret = - x[2] * x[7];
			ret += 1250000 + x[2] * x[4] - 2500;
			return ret; });

		return problem;
	}*/
	else throw std::runtime_error("Problem not implemented.\n \
			Check the documentation for the available problems.\n \
			Check if the dim template parameter is consistent with the specified problem.");
}

template <size_t dim>
double TestProblems::get_exact_value(ProblemName p)
{
	if (p == TOWNSEND && dim == 2)
	{
		return -2.0239884;
	}
	else if (p == GOMEZ_LEVY && dim == 2)
	{
		return -1.031628453;
	}
	/*else if (p == _G10 && dim == 8)
	{
		return 7049.248;
	}*/
	else
		throw std::runtime_error("Problem not implemented.\n \
			Check the documentation for the available problems.\n \
			Check if the dim template parameter is consistent with the specified problem.");
}

template <size_t dim>
RealVector<dim> TestProblems::get_exact_position(ProblemName p)
{
	if (p == TOWNSEND && dim == 2)
	{
		return RealVector<dim>(2.0052938, 1.1944509);
	}
	else if (p == GOMEZ_LEVY && dim == 2)
	{
		return RealVector<dim>(0.089842010, -0.7126564);
	}
	/*else if (p == _G10 && dim == 8)
	{
		return RealVector<dim>(579.3167, 1359.943, 5110.71, 182.0714, 295.5985, 217.979, 286.4162, 395.5979);
	}*/
	else
		throw std::runtime_error("Problem not implemented.\n \
			Check the documentation for the available problems.\n \
			Check if the dim template parameter is consistent with the specified problem.");
}
