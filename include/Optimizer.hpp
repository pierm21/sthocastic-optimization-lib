#pragma once

#include <iostream>

#include "TypeTraits.hpp"
#include "Problem.hpp"

/**
 * @brief Virtual class defining the interface for an optimization algorithm that can be employes for constrained or unconstrained optimization problems.
 *
 * @tparam dim the dimension of the space in which the objective function is defined
 */
template <size_t dim>
class Optimizer
{
protected:
	Problem<dim> problem_;

public:
	Optimizer() = delete;

	/**
	 * @brief Construct a new Optimizer object for the given problem.
	 *
	 * @param problem the Problem to be optimized
	 */
	Optimizer(const Problem<dim> &problem)
		: problem_(problem) {}
	virtual ~Optimizer() = default;

	/**
	 * @brief Initialize the optimizator to start the optimization process
	 *
	 */
	virtual void initialize() = 0;

	/**
	 * @brief Initialize the optimizator to start the optimization process exploiting parallelism
	 * 
	 */
	virtual void initialize_parallel() = 0;

	/**
	 * @brief Optimize the objective function
	 *
	 */
	virtual void optimize() = 0;

	/**
	 * @brief Optimize the objective function exploiting parallelism
	 *
	 */
	virtual void optimize_parallel() = 0;


	/**
	 * @brief Print the results of the optimization process
	 *
	 * @param out the output stream where the results will be printed
	 */
	virtual void print_results(std::ostream &out = std::cout) = 0;

	/**
	 * @brief Get the global best value object
	 *
	 * @return double the global best value
	 */
	virtual double get_global_best_value() = 0;

	/**
	 * @brief Get the global best position object
	 *
	 * @return const RealVector<dim>& a const reference to the global best position vector
	 */
	virtual const RealVector<dim> &get_global_best_position() = 0;

	/**
	 * @brief Check if the current solution is feasible according to the problem constraints (if provided)
	 *
	 * @return true if the solution is feasible
	 * @return false otherwise
	 */
	virtual bool is_feasible_solution() = 0;
};