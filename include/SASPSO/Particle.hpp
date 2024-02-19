#pragma once

#include <iostream>
#include <array>
#include <random>
#include <memory>

#include "TypeTraits.hpp"

using namespace type_traits;

/**
 * @brief Class that represents a single particle in the swarm
 *
 * @param problem_ a shared pointer to the problem to be optimized
 * @param random_generator_ a shared pointer to the marsenne twister random number generator
 * @param position_ the position of the particle in the search space
 * @param velocity_ the velocity of the particle in the search space
 * @param best_position_ the best position ever reached by the particle
 * @param best_value_ the value of the fitness function in the global best position
 * @param best_constraint_violation_ the total constraint violation in the global best position
 * @param omega_s_ the starting value of the inertia weight
 * @param omega_f_ the final value of the inertia weight
 * @param phi1_s_ the starting value of the cognitive parameter
 * @param phi1_f_ the final value of the cognitive parameter
 * @param phi2_s_ the starting value of the social parameter
 * @param phi2_f_ the final value of the social parameter
 *
 * @tparam dim the dimension of the space in which the function is defined
 */
template <size_t dim>
class Particle
{
private:
	std::shared_ptr<Problem<dim>> problem_;
	std::shared_ptr<std::mt19937> random_generator_;

	RealVector<dim> position_;
	RealVector<dim> velocity_;
	double constraint_violation_;

	RealVector<dim> best_position_;
	double best_value_;
	double best_constraint_violation_;

	double omega_s_, omega_f_;
	double phi1_s_, phi1_f_;
	double phi2_s_, phi2_f_;

public:
	/**
	 * @brief Construct a new Particle object
	 *
	 * @param problem shared pointer to the problem to be optimized
	 */
	Particle(const std::shared_ptr<Problem<dim>> &problem,
			 const std::shared_ptr<std::mt19937> &random_generator,
			 const double omega_s, const double omega_f,
			 const double phi1_s, const double phi1_f,
			 const double phi2_s, const double phi2_f)
		: problem_(problem),
		  random_generator_(random_generator),
		  omega_s_(omega_s), omega_f_(omega_f),
		  phi1_s_(phi1_s), phi1_f_(phi1_f),
		  phi2_s_(phi2_s), phi2_f_(phi2_f){};

	Particle() = default;
	~Particle() = default;

	/**
	 * @brief Initialize the particle parameters
	 */
	void initialize();

	/**
	 * @brief Update velocity, position and total constraint violation of the particle at current iteration
	 */
	void update(const RealVector<dim> &global_best_position, double violation_threshold, int iteration, int max_iter, double tol=1e-8);

	/**
	 * @brief Print the particle parameters and actual state
	 *
	 * @param out the output stream where to write the data
	 */
	void print(std::ostream &out = std::cout) const;

	/**
	 * @brief Get the actual position of the particle
	 *
	 * @return const RealVector& the position vector of the particle
	 */
	const RealVector<dim> &get_position() const { return position_; }
	/**
	 * @brief Get the best position ever reached by the particle
	 *
	 * @return const RealVector& the best position vector of the particle
	 */
	const RealVector<dim> &get_best_position() const { return best_position_; }
	/**
	 * @brief Get the best fitness value ever reached by the particle.
	 * This is the value of the fitness function in the best position.
	 *
	 * @return double the fitness value of the best position
	 */
	double get_best_value() const { return best_value_; }

	/**
	 * @brief Get the current total constraint violation of the particle
	 *
	 * @return double total constraint violation
	 */
	double get_constraint_violation() const { return constraint_violation_; }

	/**
	 * @brief Get the constraint violation in the best position found by the particle
	 *
	 * @return double best total constraint violation
	 */
	double get_best_constraint_violation() const { return best_constraint_violation_; }

	/**
	 * @brief Check if the particle is better according to the feasibility-based rule
	 *
	 * @param other the particle to be compared with
	 * @param violation_threshold the total contraint violation threshold to be used in the comparison
	 * @param tol the tolerance to be used in the comparison
	 * @return true if the particle is better than the other
	 * @return false otherwise
	 */
	bool is_better_than(const Particle<dim> &other, double violation_threshold, double tol) const;

private:
	/**
	 * @brief Utility to update the total constraint violaton at the current position
	 */
	void update_constraint_violation();

	/**
	 * @brief Utility to get the best position between the given two, following the feasibility-based rule
	 *
	 * @param value1 the fitness value on first position
	 * @param value2 the fitness value on second position
	 * @param viol1 the constraint violation on the first position
	 * @param viol2 the constraint violation on the second position
	 * @param violation_threshold the total contraint violation threshold to be used in the comparison
	 * @param tol the tolerance to be used in the comparison
	 * @return true if the first position is better than the second one
	 * @return false otherwise
	 */
	bool feasibility_rule(double value1, double value2, double viol1, double viol2, double violation_threshold, double tol) const;
};

#include "SASPSO/Particle.cpp"