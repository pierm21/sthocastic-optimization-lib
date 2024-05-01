#pragma once

#include <memory>

#include "Optimizer.hpp"
#include "Problem.hpp"
#include "SASPSO/SASPSO.hpp"
#include "ABC/ABC.hpp"

/**
 * @brief Factory class to create a specific optimizer for the given problem
 *
 * @tparam dim the dimension of the space in which the problem is defined
 */
template <std::size_t dim>
class OptimizerFactory
{
public:
	/**
	 * @brief Enumerator to represent the available optimizers
	 * - SelfAdaptiveSPSO: Self Adaptive Particle Swarm Optimization
	 */
	enum OptimizerName
	{
		SelfAdaptiveSPSO,
		ArtificialBeeColony
	};

	/**
	 * @brief Method to create a specific optimizer for the given problem given the optimizer name.
	 * Throws runtime error if the name is not available.
	 * The optimizer is created with its default parameters.
	 *
	 * @param name the name of the choosen optimizer
	 * @param problem the problem to be optimized
	 * @return std::unique_ptr<Optimizer<dim>> the pointer to the created optimizer
	 */
	std::unique_ptr<Optimizer<dim>> create(OptimizerName name, Problem<dim>& problem)
	{
		if(name == OptimizerName::SelfAdaptiveSPSO)
			return std::make_unique<SASPSO<dim>>(problem);
		else if(name == OptimizerName::ArtificialBeeColony)
			return std::make_unique<ABC<dim>>(problem);
		else
			throw std::runtime_error("Invalid optimizer name");
	}
};
