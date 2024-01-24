#pragma once

#include "SASPSO/SASPSO.hpp"

template <std::size_t dim>
void SASPSO<dim>::initialize()
{
	// Initialize the global best variables
	global_best_value_ = std::numeric_limits<double>::max();
	global_best_position_ = RealVector<dim>::Zero();
	global_best_constraint_violation_ = std::numeric_limits<double>::max();

	// Instantiate the marsenne twister
	std::random_device rand_dev;
	std::shared_ptr<std::mt19937> generator = std::make_shared(rand_dev());
	// Create a shared pointer to the problem
	std::shared_ptr<Problem<dim>> problem = std::make_shared(problem_);

	// Swarm creation
	for(size_t i = 0; i < swarm_size_; ++i)
	{
		// Create the particle and initialize it
		swarm_.push_back(Particle<dim>(problem, generator, omega_s_, omega_f_, phi1_s_, phi1_f_, phi2_s_, phi2_f_));
		swarm_[i].initialize();
		// Update the global best variables
		if(swarm_[i].get_best_value() < global_best_value_)
		{
			//global_best_value_ = swarm_[i].get_best_value();
			//global_best_position_ = swarm_[i].get_best_position();
			//global_best_constraint_violation_ = swarm_[i].get_best_constraint_violation();
		}
	}
}

template <std::size_t dim>
void SASPSO<dim>::optimize()
{
	std::cout << "SASPSO::optimize()" << std::endl; //TODO: implement
}

template <std::size_t dim>
void SASPSO<dim>::print_results(std::ostream &out)
{
	std::cout << "SASPSO::print_results()" << std::endl; //TODO: implement
}

template <std::size_t dim>
double SASPSO<dim>::get_global_best_value()
{
	return 0.0; //TODO: implement
}

template <std::size_t dim>
const RealVector<dim> &SASPSO<dim>::get_global_best_position()
{
	return RealVector<dim>::Zero(dim,1); //TODO: implement
}

template <std::size_t dim>
bool SASPSO<dim>::is_feasible_solution()
{
	return false; //TODO: implement
}