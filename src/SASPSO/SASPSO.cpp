#pragma once

#include "SASPSO/SASPSO.hpp"

template <std::size_t dim>
void SASPSO<dim>::initialize()
{
	// Instantiate the marsenne twister
	std::random_device rand_dev;
	std::shared_ptr<std::mt19937> generator = std::make_shared<std::mt19937>(rand_dev());
	// Create a shared pointer to the problem
	std::shared_ptr<Problem<dim>> problem = std::make_shared<Problem<dim>>(Optimizer<dim>::problem_);

	// Initialize the array of total constraint violations
	std::vector<double> total_violations;

	// Initialize the first particle
	swarm_.emplace_back(problem, generator, omega_s_, omega_f_, phi1_s_, phi1_f_, phi2_s_, phi2_f_);
	swarm_[0].initialize();

	// Initialize the global best
	global_best_index_ = 0;

	// Initialize the remaining particle in the swarm
	for (std::size_t i = 1; i < swarm_size_; ++i)
	{
		// Create and initialize the particle
		swarm_.emplace_back(problem, generator, omega_s_, omega_f_, phi1_s_, phi1_f_, phi2_s_, phi2_f_);
		swarm_[i].initialize();
		// Update the global best
		if (swarm_[i].is_better_than(swarm_[global_best_index_], tol_))
			global_best_index_ = i;
		// Add the constraint violation to the array
		total_violations.push_back(swarm_[i].get_best_constraint_violation());
	}

	// Initialize the violation threshold as the median of the total violations
	std::sort(total_violations.begin(), total_violations.end());
	if (swarm_size_ % 2 == 0)
		violation_threshold_ = (total_violations[swarm_size_ / 2] + total_violations[swarm_size_ / 2 - 1]) / 2;
	else
		violation_threshold_ = total_violations[swarm_size_ / 2];
}

template <std::size_t dim>
void SASPSO<dim>::optimize()
{
	int current_iter = 0;
    double temp_value = 0.0;
	int feasible_particles = 0;

    // Outer optimization loop over all the iterations
	while (current_iter < max_iter_)
    {
		// Reset the number of feasible particles for the current iteration
		feasible_particles = 0;

		// Process each particle of the swarm
        for (size_t i = 0; i < swarm_size_; ++i)
		{
			// Update the particle
		   	swarm_[i].update(swarm_[global_best_index_].get_best_position(), current_iter, max_iter_, tol_);
			// Check if the particle is feasible
			if (swarm_[i].get_best_constraint_violation() <= violation_threshold_) {
				feasible_particles++;
				// Update global best position
				if (swarm_[i].is_better_than(swarm_[global_best_index_], tol_))
					global_best_index_ = i;
			}
        }

		// Update the violation threshold according to the proportion of feasible particles
		violation_threshold_ = violation_threshold_ * (1 - (feasible_particles / (double)swarm_size_));

		std::cout << current_iter << " | " << swarm_[global_best_index_].get_best_value() << " | " << swarm_[global_best_index_].get_best_constraint_violation() << " | " << feasible_particles << " | " << violation_threshold_ << std::endl;

		// Update the current iteration
        current_iter++;
    }
}

template <std::size_t dim>
void SASPSO<dim>::print_results(std::ostream &out)
{
	if(global_best_index_ == -1)
	{
		out << "No results to show" << std::endl;
		return;
	}
	out << "Best value: " << swarm_[global_best_index_].get_best_value() << std::endl;
	out << "Best position: " << swarm_[global_best_index_].get_best_position() << std::endl;
	out << "Total constraint violation: " << swarm_[global_best_index_].get_best_constraint_violation() << std::endl;
}

template <std::size_t dim>
double SASPSO<dim>::get_global_best_value()
{
	return swarm_[global_best_index_].get_best_value();
}

template <std::size_t dim>
const RealVector<dim> &SASPSO<dim>::get_global_best_position()
{
	return swarm_[global_best_index_].get_best_position();
}

template <std::size_t dim>
bool SASPSO<dim>::is_feasible_solution()
{
	return swarm_[global_best_index_].get_best_constraint_violation() <= tol_;
}