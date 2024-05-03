#pragma once

#include <iostream>
#include <array>
#include <random>
#include <memory>

#include "TypeTraits.hpp"
#include "Optimizer.hpp"
#include "Particle.hpp"


using namespace type_traits;

/**
 * @brief Class that represents a single Bee in the colony
 *
 * @param fi_value_ the value of the fitness function in the global best position
 * @param fi_probability_ proportional to the fitness value of the Bee, used in updating the position
 * @param failure_counter_ the number of iterations in which the Bee has not improved
 * @param constraint_violation_ the sontraint violation of the Bee's position
 *
 * @tparam dim the dimension of the space in which the function is defined
 */
template <size_t dim>
class Bee : public Particle<dim>
{
private:
	double value_;
	double fi_probability_ ;
	int failure_counter_ ;
	double constraint_violation_;

public:
	/**
	 * @brief Construct a new Bee object
	 *
	 * @param problem shared pointer to the problem to be optimized
	 */
	Bee(const std::shared_ptr<Problem<dim>> &problem,
			 const std::shared_ptr<std::mt19937> &random_generator)
		: Particle<dim> (problem, random_generator){};

	Bee() = default;
	~Bee() = default;

	/**
	 * @brief Initialize the Bee parameters
	 */
	void initialize() override;

	/**
	 * @brief Update the Bee position, using the three different strategies: employer, onlooker and scout
	 * 
	 */
    void update(const int limit, const double MR, double violation_threshold, const std::vector<Bee<dim>>& colony_, double tol=1e-8);

	/**
	 * @brief Print the Bee parameters and actual state
	 *
	 * @param out the output stream where to write the data
	 */
	void print(std::ostream &out = std::cout) const override;

	/**
	 * @brief Get the fitness value for the acutal position occupied by the Bee.
	 *
	 * @return double the fitness value of the best position
	 */
	double get_value() const { return value_; }

	/**
	 * @brief Get the failure counter of the Bee
	 *
	 * @return int the failure counter
	 */
	int get_failure_counter() const { return failure_counter_; }

	/**
	 * @brief Get the current total constraint violation of the particle
	 *
	 * @return double total constraint violation
	 */
	double get_constraint_violation() const { return constraint_violation_; }


private:
	/**
	 * @brief Utility to update the total constraint violaton at the current position
	 * 
	 * @return double indicating the amount of constraint violation 
	 */
	double compute_constraint_violation(const RealVector<dim> &position) const;

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

#include "ABC/Bee.cpp"
