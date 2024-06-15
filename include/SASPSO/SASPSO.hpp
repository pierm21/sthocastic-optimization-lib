#pragma once

#include "Optimizer.hpp"
#include "Particle_SASPSO.hpp"

/**
 * @brief Specialization of the Optimizer class implementing the Self Adaptive Standard Particle Swarm Optimization (SASPSO) algorithm.
 * All the information regarding the analysis and the theory behind this algorithm can be found in the following paper: http://dx.doi.org/10.1155/2016/8627083
 *
 * @param swarm_size the number of particles in the swarm
 * @param max_iter the maximum number of iterations
 * @param tol the tolerance used for checking constraint conditions
 * @param omega_s the starting value of the inertia weight
 * @param omega_f the final value of the inertia weight
 * @param phi1_s the starting value of the cognitive parameter
 * @param phi1_f the final value of the cognitive parameter
 * @param phi2_s the starting value of the social parameter
 * @param phi2_f the final value of the social parameter
 * @param swarm_ the std::vector storing the swarm of Particle_SASPSO instances
 * @param global_best_index_ the index in the swarm array of the global best particle
 * @param violation_threshold the threshold used for considering a particle as feasible or not
 *
 * @tparam dim the dimension of the space in which the function is defined
 */
template <std::size_t dim>
class SASPSO : public Optimizer<dim>
{
private:
	int swarm_size_;
	int max_iter_;
	double tol_;

	double omega_s_, omega_f_;
	double phi1_s_, phi1_f_;
	double phi2_s_, phi2_f_;

	std::vector<Particle_SASPSO<dim>> swarm_;

	size_t global_best_index_;

	double violation_threshold_;

public:
	/**
	 * @brief Construct a new SASPSO optimizator object for the given problem.
	 *
	 * @param problem the Problem to be optimized
	 * @param swarm_size the number of particles in the swarm
	 * @param max_iter the maximum number of iterations to be performed
	 * @param tol the tolerance used for checking constraint conditions
	 * @param omega_s the starting value of the inertia weight
	 * @param omega_f the final value of the inertia weight
	 * @param phi1_s the starting value of the cognitive parameter
	 * @param phi1_f the final value of the cognitive parameter
	 * @param phi2_s the starting value of the social parameter
	 * @param phi2_f the final value of the social parameter
	 */
	SASPSO(const Problem<dim> &problem,
		   int swarm_size = 100, int max_iter = 2000, double tol = 1e-6,
		   const double omega_s = 0.9, const double omega_f = 0.4,
		   const double phi1_s = 2.5, const double phi1_f = 0.3,
		   const double phi2_s = 0.3, const double phi2_f = 2.5)
		: Optimizer<dim>(problem),
		  swarm_size_(swarm_size), max_iter_(max_iter), tol_(tol),
		  omega_s_(omega_s), omega_f_(omega_f),
		  phi1_s_(phi1_s), phi1_f_(phi1_f),
		  phi2_s_(phi2_s), phi2_f_(phi2_f),
		  global_best_index_(-1){};

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
	 * @brief Optimize the given problem and store the history of the best value and violation found every interval iterations.
	 * It can also log the time evolution of all the particles to a file, if provided.
	 *
	 * @param optimum_history the vector where to store the history of the best value found
	 * @param violation_history the vector where to store the history of the best contraint violation found
	 * @param out a pointer to the output stream where to print the time evolution of all the particles (default: nullptr)
	 * @param interval the sampling interval in number of iterations
	 */
	void optimize(std::vector<double> &optimum_history, std::vector<double> &violation_history, const int interval = 50, std::ostream *out = nullptr);

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
	void print_results(std::ostream &out = std::cout) const override;

	/**
	 * @brief Get the actual global best value found by the algorithm
	 *
	 * @return double the global best value
	 */
	double get_global_best_value() const override;

	/**
	 * @brief Get the position of the global best minimum found by the algorithm
	 *
	 * @return const RealVector<dim>& a const reference to the global best position vector
	 */
	const RealVector<dim> &get_global_best_position() const override;

	/**
	 * @brief Checks if the global best found by the algorithm is a feasible solution according to the constraints
	 *
	 * @return true if the global best is a feasible solution
	 * @return false otherwise
	 */
	bool is_feasible_solution() const;
};

#include "SASPSO/SASPSO.cpp"
