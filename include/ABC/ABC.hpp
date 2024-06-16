#pragma once

#include "Optimizer.hpp"
#include "Bee.hpp"
#include "mpi.h"

/**
 * @brief Specialization of the Optimizer class implementing the Artificial Bee Colony Optimization (ABC) algorithm.
 * All the information regarding the analysis and the theory behind this algorithm can be found in the following paper:
 * https://link.springer.com/content/pdf/10.1007/978-3-540-72950-1_77.pdf
 *
 * @param colony_size the number of particles in the swarm
 * @param max_iter the maximum number of iterations
 * @param tol the tolerance used for checking constraint conditions
 * @tparam dim the dimension of the space in which the function is defined
 */
template <std::size_t dim>
class ABC : public Optimizer<dim>
{
private:
	int colony_size_;
	int max_iter_;
	int limit_;
	double MR_;
	int SPP_;

	RealVector<dim> global_best_position_;
	double global_best_value_;
	double global_best_constraint_violation_;

	std::vector<Bee<dim>> colony_;

	double violation_threshold_;
	double tol_;

public:
	/**
	 * @brief Construct a new ABC optimizator object for the given problem.
	 *
	 * @param problem the Problem to be optimized
	 * @param colony_size the number of particles in the swarm
	 */

	ABC(const Problem<dim> &problem, int colony_size = 40, int max_iter = 6000, int limit = -1, int SPP = -1, double MR = 0.8, double tol = 1e-6)
		: Optimizer<dim>(problem),
		  colony_size_(colony_size), max_iter_(max_iter),
		  MR_(MR), tol_(tol),
		  limit_(limit == -1 ? static_cast<int>(colony_size * dim * 0.5) : limit), // TODO: esnure that dim can be used here
		  SPP_(SPP == -1 ? static_cast<int>(colony_size * dim * 0.5) : SPP){};

	/**
	 * @brief Initialize the optimizator to start the optimization process
	 * This method initializes the swarm of particles and all the parameters needed by the algorithm
	 */
	void initialize() override;

	/**
	 * @brief Initialize the optimizator to start the optimization process using OMP parallel constructs
	 * This initialization must be used if the optimization will be performed using optimize_parallel method
	 */
	void initialize_parallel() override;

	/**
	 * @brief Optimize the given problem
	 */
	void optimize() override;

	/**
	 * @brief Optimize the given problem and store the history of the best value found every interval iterations
	 *
	 * @param optimum_history the vector where to store the history of the best value found
	 * @param violation_history the vector where to store the history of the best contraint violation found
	 * @param interval the sampling interval in number of iterations
	 */
	void optimize(std::vector<double> &optimum_history, std::vector<double> &violation_history, std::vector<double> &feasible_history, const int interval = 50 , std::ostream *out = nullptr) override;

	/**
	 * @brief Optimize the given problem using OMP thread level parallel constructs
	 */
	void optimize_parallel() override;

	/**
	 * @brief Optimize the given problem using OMP thread level parallel constructs and store the history of the best value found every interval iterations
	 *
	 * @param optimum_history the vector where to store the history of the best value found
	 * @param violation_history the vector where to store the history of the best contraint violation found
	 * @param feasible_history the vector where to store the history of the number of feasible solutions found
	 * @param interval the sampling interval in number of iterations
	 */
	void optimize_parallel(std::vector<double> &optimum_history, std::vector<double> &violation_history, std::vector<double> &feasible_history, const int interval = 50) override;

	/**
	 * @brief Print the results of the optimization process to the given output stream
	 *
	 * @param out output stream where to print results
	 */
	void print_results(std::ostream &out = std::cout) const override;

	/**
	 * @brief Get the actual global best value found by the algorithm
	 *
	 * @return double the global best value
	 */
	double get_global_best_value() const override { return global_best_value_; }

	/**
	 * @brief Get the position of the global best minimum found by the algorithm
	 *
	 * @return const RealVector<dim>& a const reference to the global best position vector
	 */
	const RealVector<dim> &get_global_best_position() const override { return global_best_position_; }

	double get_global_best_constraint_violation() const {return global_best_constraint_violation_;};
};

#include "ABC/ABC.cpp"
