#pragma once

#include "Optimizer.hpp"
#include "Particle.hpp"

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
 * @param swarm_ the std::vector storing the swarm of Particle instances
 * @param global_best_constraint_violation_ the constraint violation of the global best particle
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

	std::vector<Particle<dim>> swarm_;

	double global_best_constraint_violation_;
	double violation_threshold;

public:
	/**
	 * @brief Construct a new SASPSO optimizator object for the given problem.
	 *
	 * @param problem the Problem to be optimized
	 * @param swarm_size the number of particles in the swarm
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
		   const double omega_s = 0.9, const double omega_f = 0.1,
		   const double phi1_s = 2.5, const double phi1_f = 0.1,
		   const double phi2_s = 0.1, const double phi2_f = 2.5)
		: Optimizer<dim>(problem),
		  swarm_size_(swarm_size), max_iter_(max_iter), tol_(tol),
		  omega_s_(omega_s), omega_f_(omega_f),
		  phi1_s_(phi1_s), phi1_f_(phi1_f),
		  phi2_s_(phi2_s), phi2_f_(phi2_f){};

	/**
	 * @brief Initialize the optimizator to start the optimization process
	 * This method initializes the swarm of particles and all the parameters needed by the algorithm.
	 */
	void initialize() override;
	/**
	 * @brief Optimize the given problem
	 */
	void optimize() override;
	/**
	 * @brief Print the results of the optimization process to the given output stream
	 *
	 * @param out output stream where to print results
	 */
	void print_results(std::ostream &out = std::cout) override;
};

#include "SASPSO/SASPSO.hpp"
