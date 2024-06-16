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
	static std::unique_ptr<Optimizer<dim>> create(
		OptimizerName name, Problem<dim>& problem, int size, int iter, double tol=1e-6, double omega_s = 0.9,
		 double omega_f = 0.4, double phi1_s = 2.5, double phi1_f = 0.3, double phi2_s = 2.5, double phi2_f = 0.3)
	{
		if(name == OptimizerName::SelfAdaptiveSPSO)
			return std::make_unique<SASPSO<dim>>(problem, size, iter, tol);
		else if(name == OptimizerName::ArtificialBeeColony)
			return std::make_unique<ABC<dim>>(problem, size, iter);
		else
			throw std::runtime_error("Invalid optimizer name");
	}

	/**
	 * @brief Get the string name object
	 *
	 * @param name the enum type of the optimizer
	 * @return std::string the string name of the optimizer
	 */
	static const std::string get_string_name(const OptimizerName name)
	{
		if(name == OptimizerName::SelfAdaptiveSPSO)
			return "saspso";
		else if(name == OptimizerName::ArtificialBeeColony)
			return "abc";
		else
			throw std::runtime_error("Invalid optimizer name");
	}
};
