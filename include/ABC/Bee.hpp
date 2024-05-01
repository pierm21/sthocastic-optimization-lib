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
 * @param violation_ the sontraint violation of the Bee's position
 *
 * @tparam dim the dimension of the space in which the function is defined
 */
template <size_t dim>
class Bee : public Particle<dim>
{
private:
	double fi_value_;
	double fi_probability_ ;
	int failure_counter_ ;
	double violation_;

public:
	/**
	 * @brief Construct a new Bee object
	 *
	 * @param problem shared pointer to the problem to be optimized
	 */
	Bee(const std::shared_ptr<Problem<dim>> &problem,
			 const std::shared_ptr<std::mt19937> &random_generator)
		: Particle<dim> (problem, random_generator),
		  failure_counter_(0){};

	Bee() = default;
	~Bee() = default;

	/**
	 * @brief Initialize the Bee parameters
	 */
	void initialize();

	/**
	 * @brief Update the Bee position according to the empoyer exploitation strategy
	 */
	void employer_update();

	/**
	 * @brief Update the Bee position according to the onlooker exploitation strategy
	 */
	void outlooker_update();

	/**
	 * @brief Update the Bee position according to the scout exploration strategy
	 */
	void scout_update();

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
	double get_value() const { return fi_value_; }

	/**
	 * @brief Get the failure counter of the Bee
	 *
	 * @return int the failure counter
	 */
	int get_failure_counter() const { return failure_counter_; }

};

//#include "SASPSO/Bee.cpp"
