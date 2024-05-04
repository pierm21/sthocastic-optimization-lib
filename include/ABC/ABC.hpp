#pragma once

#include "Optimizer.hpp"
#include "Bee.hpp"

/**
 * @brief Specialization of the Optimizer class implementing the Artificial Bee Colony Optimization (ABC) algorithm.
 * All the information regarding the analysis and the theory behind this algorithm can be found in the following paper:
 * https://link.springer.com/content/pdf/10.1007/978-3-540-72950-1_77.pdf
 *
 * @param colony_size the number of particles in the swarm
 * @param max_iter the maximum number of iterations
 * @param tol the tolerance used for checking constraint conditions
 * @param global_best_index_ the index in the swarm array of the global best particle
 *
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

	int global_best_index_;

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
	
	ABC(const Problem<dim> &problem, int limit, int SPP, int colony_size = 40, int max_iter = 6000, double MR = 0.8, double tol = 1e-6)
		: Optimizer<dim>(problem),
		  colony_size_(colony_size), max_iter_(max_iter), 
		  global_best_index_(-1),MR_(MR),
		  limit_(limit), SPP_(SPP), tol_(tol){};

	ABC(const Problem<dim> &problem, int colony_size = 40, int max_iter = 6000, double MR = 0.8, double tol = 1e-6)
		: Optimizer<dim>(problem),
		  colony_size_(colony_size), max_iter_(max_iter), 
		  global_best_index_(-1),MR_(MR),
		  limit_(static_cast<int>(colony_size*dim*0.5)),				//TODO: esnure that dim can be used here
		  SPP_(static_cast<int>(colony_size*dim*0.5)){};


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
	void optimize(std::vector<double> &optimum_history, std::vector<double> &violation_history, const int interval = 50);

	/**
	 * @brief Optimize the given problem using OMP thread level parallel constructs
	 */
	void optimize_parallel() override;

	/**
	 * @brief Optimize the given problem using OMP thread parallelism and store the history of the best value found every interval iterations
	 *
	 * @param optimum_history the vector where to store the history of the best value found
	 * @param violation_history the vector where to store the history of the best contraint violation found
	 * @param feasible_history the vector where to store the history of the number of feasible solutions found
	 * @param threshold_history the vector where to store the history of the violation thresholds
	 * @param interval the sampling interval in number of iterations
	 */
	void optimize_parallel(std::vector<double> &optimum_history, std::vector<double> &violation_history, std::vector<double> &feasible_history, std::vector<double> &threshold_history, const int interval = 50);

	/**
	 * @brief Print the results of the optimization process to the given output stream
	 *
	 * @param out output stream where to print results
	 */
	void print_results(std::ostream &out = std::cout) override;

	/**
	 * @brief Get the actual global best value found by the algorithm
	 *
	 * @return double the global best value
	 */
	double get_global_best_value() override;

	/**
	 * @brief Get the position of the global best minimum found by the algorithm
	 *
	 * @return const RealVector<dim>& a const reference to the global best position vector
	 */
	const RealVector<dim> &get_global_best_position() override;

	/**
	 * @brief Checks if the global best found by the algorithm is a feasible solution according to the constraints
	 *
	 * @return true if the global best is a feasible solution
	 * @return false otherwise
	 */
	bool is_feasible_solution() override;
};

#include "ABC/ABC.cpp"
