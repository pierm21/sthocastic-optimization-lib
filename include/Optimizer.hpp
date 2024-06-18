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
	 * @brief Optimize the objective function
	 *
	 */
	virtual void optimize() = 0;

	/**
	 * @brief Optimize the objective fufunction and store the results in the given vectors
	 *
	 */
	virtual void optimize(std::vector<double> &optimum_history, std::vector<double> &violation_history, std::vector<double> &feasible_history, const int interval = 50, std::ostream *out = nullptr) = 0;

	/**
	 * @brief Optimize the objective function exploiting parallelism
	 *
	 */
	virtual void optimize_parallel() = 0;

	/**
	 * @brief 
	 * 
	 * @param optimum_history 
	 * @param violation_history 
	 * @param feasible_history 
	 * @param interval
	 * @param out
	 */
	virtual void optimize_parallel(std::vector<double> &optimum_history, std::vector<double> &violation_history, std::vector<double> &feasible_history, const int interval = 50) = 0;


	/**
	 * @brief Print the results of the optimization process
	 *
	 * @param out the output stream where the results will be printed
	 */
	virtual void print_results(std::ostream &out = std::cout) const = 0;

	/**
	 * @brief Get the global best value object
	 *
	 * @return double the global best value
	 */
	virtual double get_global_best_value() const = 0;

	/**
	 * @brief Get the global best position object
	 *
	 * @return const RealVector<dim>& a const reference to the global best position vector
	 */
	virtual const RealVector<dim> &get_global_best_position() const = 0;

	/**
	 * @brief Set the verbosity for logging purposes
	 *
	 * @param log_verbose the verbosity for logging purposes
	 */
	void set_log_verbose(bool log_verbose) { log_verbose_ = log_verbose; }

};