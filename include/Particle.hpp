#pragma once

#include <iostream>
#include <array>
#include <random>
#include <memory>

#include "TypeTraits.hpp"
#include "Optimizer.hpp"


using namespace type_traits;

/**
 * @brief Virtual Class that represents a single particle in the swarm/colony
 *
 * @param problem_ a shared pointer to the problem to be optimized
 * @param random_generator_ a shared pointer to the marsenne twister random number generator
 * @param position_ the position of the bee in the search space
 *
 * @tparam dim the dimension of the space in which the function is defined
 */
template <size_t dim>
class Particle
{
protected:
	std::shared_ptr<Problem<dim>> problem_;
	std::shared_ptr<std::mt19937> random_generator_;

	RealVector<dim> position_;

public:
	/**
	 * @brief Construct a new Particle object
	 *
	 * @param problem shared pointer to the problem to be optimized
	 */

	Particle(const std::shared_ptr<Problem<dim>> &problem,
			 const std::shared_ptr<std::mt19937> &random_generator)
		: problem_(problem),
		  random_generator_(random_generator){};

	virtual ~Particle() = default;
    Particle() = default;


	/**
	 * @brief Initialize the Particle parameters
	 */
	virtual void initialize() = 0;

	/**
	 * @brief Print the Particle parameters and actual state
	 *
	 * @param out the output stream where to write the data
	 */
	virtual void print(std::ostream &out = std::cout) const = 0;

	/**
	 * @brief Get the actual position of the Particle
	 *
	 * @return const RealVector& the position vector of the Particle
	 */
	const RealVector<dim> &get_position() const { return position_; }

};
