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
	bool log_verbose_;

public:
	Optimizer() = delete;

	/**
	 * @brief Construct a new Optimizer object for the given problem.
	 *
	 * @param problem the Problem to be optimized
	 * @param log_verbose a boolean to enable/disable the verbosity for logging purposes
	 */
	Optimizer(const Problem<dim> &problem, bool log_verbose = false)
		: problem_(problem), log_verbose_(log_verbose) {}
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
	 * @brief Optimize the objective function and print to the given streams the history of the optimization process if log_verbose is true
	 *
	 * @param optimum_history the stream to which print the history as: iteration, fitness, constraint_violation, feasible_particles
	 * @param simulation_history the stream to which print the data to produce a simulation of the optimization process
	 * @param interval the interval in number of iterations to print the history
	 */
	virtual void optimize(std::ostream& history_out = std::cout, std::ostream& simulation_out = std::cout, const int interval = 50) = 0;

	/**
	 * @brief Optimize the objective function in parallel and print to the given streams the history of the optimization process if log_verbose is true
	 * The parallelism type depend on the specific implementation of the optimizer. May be multi-threading only or multi-processing + multi-threading
	 *
	 * @param optimum_history the stream to which print the history as: iteration, fitness, constraint_violation, feasible_particles
	 * @param simulation_history the stream to which print the data to produce a simulation of the optimization process
	 * @param interval the interval in number of iterations to print the history
	 */
	virtual void optimize_parallel(std::ostream& history_out = std::cout, std::ostream& simulation_out = std::cout, const int interval = 50) = 0;

	/**
	 * @brief Print the results of the optimization process
	 *
	 * @param out the output stream where the results will be printed
	 */
	virtual void print_results(std::ostream &out = std::cout) const = 0;

	/**
	 * @brief Get the global best value object
	 *
	 * @return double the global best value found up to the moment
	 */
	virtual double get_global_best_value() const = 0;

	/**
	 * @brief Get the global best position object
	 *
	 * @return const RealVector<dim>& a const reference to the global best position vector found up to the moment
	 */
	virtual const RealVector<dim> &get_global_best_position() const = 0;

	/**
	 * @brief Set the verbosity for logging purposes
	 *
	 * @param log_verbose the verbosity for logging purposes
	 */
	void set_log_verbose(bool log_verbose) { log_verbose_ = log_verbose; }

};