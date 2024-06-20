#pragma once

#include "SASPSO/SASPSO.hpp"
#include <omp.h>
#include <iomanip>

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
	// Add the constraint violation to the array
	total_violations.push_back(swarm_[0].get_best_constraint_violation());

	// Initialize the global best
	global_best_index_ = 0;

	// Initialize the remaining particle in the swarm
	for (std::size_t i = 1; i < swarm_size_; ++i)
	{
		// Create and initialize the particle
		swarm_.emplace_back(problem, generator, omega_s_, omega_f_, phi1_s_, phi1_f_, phi2_s_, phi2_f_);
		swarm_[i].initialize();
		// Add the constraint violation to the array
		total_violations.push_back(swarm_[i].get_best_constraint_violation());
	}

	// Initialize the violation threshold as the median of the total violations
	std::sort(total_violations.begin(), total_violations.end());
	if (swarm_size_ % 2 == 0)
		violation_threshold_ = (total_violations[swarm_size_ / 2] + total_violations[swarm_size_ / 2 - 1]) / 2;
	else
		violation_threshold_ = total_violations[swarm_size_ / 2];

	// Update the global best considering the violation threshold
	for (std::size_t i = 1; i < swarm_size_; ++i)
		if (swarm_[i].is_better_than(swarm_[global_best_index_], violation_threshold_, tol_))
			global_best_index_ = i;
}

template <std::size_t dim>
void SASPSO<dim>::initialize_parallel()
{
	// Create a shared pointer to the problem
	std::shared_ptr<Problem<dim>> problem = std::make_shared<Problem<dim>>(Optimizer<dim>::problem_);

	// Initialize the array of total constraint violations. Needed since threads performs direct access to the array
	std::vector<double> total_violations;
	total_violations.resize(swarm_size_);

	// Initialize the swarm. Needed since threads performs direct access to the array
	swarm_.resize(swarm_size_);

	// Add the constraint violation to the array
	total_violations[0] = swarm_[0].get_best_constraint_violation();

	// Initialize the global best
	global_best_index_ = 0;

// Spawn the threads to initialize the particles
#pragma omp parallel shared(swarm_, total_violations, global_best_index_, problem)
	{
		// Instantiate the marsenne twister, one for each thread since it is not thread safe
		std::random_device rand_dev;
		std::shared_ptr<std::mt19937> generator = std::make_shared<std::mt19937>(rand_dev());

		// Initialize all the particles. Barrier at the end of the parallel region is needed
#pragma omp for schedule(static)
		for (std::size_t i = 0; i < swarm_size_; ++i)
		{
			// Create and initialize the particle
			swarm_[i] = Particle_SASPSO(problem, generator, omega_s_, omega_f_, phi1_s_, phi1_f_, phi2_s_, phi2_f_);
			swarm_[i].initialize();
			// Add the constraint violation to the array
			total_violations[i] = swarm_[i].get_best_constraint_violation();
		}
#pragma omp single
		{
			// Initialize the violation threshold as the median of the total violations
			std::sort(total_violations.begin(), total_violations.end());
			if (swarm_size_ % 2 == 0)
				violation_threshold_ = (total_violations[swarm_size_ / 2] + total_violations[swarm_size_ / 2 - 1]) / 2;
			else
				violation_threshold_ = total_violations[swarm_size_ / 2];
		}

		// Initialize the thread local best index
		std::size_t local_best_index = global_best_index_;

		// Find the local best index for each thread
#pragma omp for nowait schedule(static)
		for (std::size_t i = 1; i < swarm_size_; ++i)
		{
			if (swarm_[i].is_better_than(swarm_[local_best_index], violation_threshold_, tol_))
				local_best_index = i;
		}

		// Reduction: Each thread updates the global best index if it has a better solution
#pragma omp critical(update_global_best_init)
		{
			if (swarm_[local_best_index].is_better_than(swarm_[global_best_index_], violation_threshold_, tol_))
				global_best_index_ = local_best_index;
		}
	}
}

template <std::size_t dim>
void SASPSO<dim>::optimize(std::ostream& history_out, std::ostream& simulation_out, const int interval)
{
	int current_iter = 0;

	// Outer optimization loop over all the iterations
	while (current_iter < max_iter_)
	{
		// Reset the number of feasible particles for the current iteration
		int feasible_particles = 0;

		// Process each particle of the swarm
		for (size_t i = 0; i < swarm_size_; ++i)
		{
			// Update the particle
			swarm_[i].update(swarm_[global_best_index_].get_best_position(), violation_threshold_, current_iter, max_iter_, tol_);
			// Check if the particle is feasible
			if (swarm_[i].get_best_constraint_violation() <= violation_threshold_)
				feasible_particles++;
		}

		// Update the violation threshold according to the proportion of feasible particles
		violation_threshold_ = violation_threshold_ * (1.0 - (feasible_particles / (double)swarm_size_));
		violation_threshold_ = violation_threshold_ < tol_ ? 0 : violation_threshold_;

		// Update the global best position
		for (std::size_t i = 0; i < swarm_size_; ++i)
			if (swarm_[i].is_better_than(swarm_[global_best_index_], violation_threshold_, tol_))
				global_best_index_ = i;

		// Log the optimization process in both history and simulation streams.
		// NB: Branch prediction will optimize well the following branch since it has a constant condition
		if (Optimizer<dim>::log_verbose_ && current_iter % interval == 0)
		{
			// Print Iters,value,violation,threshold,feasible_particles to the history
			history_out << current_iter << ","
						<< swarm_[global_best_index_].get_best_value() << ","
						<< swarm_[global_best_index_].get_best_constraint_violation() << ","
						<< violation_threshold_ << ","
						<< feasible_particles <<  std::endl;
			// Print all the particles position and flag the best one to the simulation stream
			for (int i = 0; i < swarm_size_; i++)
			{
				simulation_out << current_iter << ",";
				for (int d = 0; d < dim; d++)
					simulation_out << swarm_[i].get_best_position()[d] << ",";
				simulation_out << ((i == global_best_index_) ? 1 : 0) << std::endl;
			}
		}

		// Update the current iteration
		current_iter++;
	}
}

template <std::size_t dim>
void SASPSO<dim>::optimize_parallel(std::ostream& history_out, std::ostream& simulation_out, const int interval)
{
	int current_iter = 0;
	std::cout << "iter" << " | " << "global best" << " | " << "global violation" << " | " << "feasible particles" << " | " << "violation threshold" << " | " << "global best index" << std::endl;

	// Outer optimization loop over all the iterations. This loop isn't parallelized, the inner loop is
	// Even if the thread creation is inside the loop, the compiler optimizes it and creates the threads only once (same execution time)
	while (current_iter < max_iter_)
	{
		// Reset the number of feasible particles for the current iteration
		int feasible_particles = 0;

		// Initialize the threads
#pragma omp parallel shared(global_best_index_, swarm_, feasible_particles, violation_threshold_)
		{
			// Initialize the thread local best index
			int local_best_index = global_best_index_;

			// Process each particle subdividing them among all threads
#pragma omp for schedule(static) reduction(+ : feasible_particles)
			for (size_t i = 0; i < swarm_size_; ++i)
			{
				//  Update the particle
				swarm_[i].update(swarm_[global_best_index_].get_best_position(), violation_threshold_, current_iter, max_iter_, tol_);
				// Check if the particle is feasible
				if (swarm_[i].get_best_constraint_violation() <= violation_threshold_)
					feasible_particles++;
			} // implicit barrier

			// Update the violation threshold according to the proportion of feasible particles
#pragma omp single
			{
				violation_threshold_ = violation_threshold_ * (1 - (feasible_particles / (double)swarm_size_));
				violation_threshold_ = violation_threshold_ < tol_ ? 0 : violation_threshold_;
			} // implicit barrier

			// Update local best position
#pragma omp for nowait schedule(static)
			for (size_t i = 0; i < swarm_size_; ++i)
			{
				if (swarm_[i].is_better_than(swarm_[local_best_index], violation_threshold_, tol_))
					local_best_index = i;
			}

			// Each thread updates the global best index if it has a better solution
#pragma omp critical(update_global_best_opt)
			{
				if (swarm_[local_best_index].is_better_than(swarm_[global_best_index_], violation_threshold_, tol_))
					global_best_index_ = local_best_index;
			} // implicit barrier
		}

		// Log the optimization process in both history and simulation streams.
		// NB: Branch prediction will optimize well the following branch since it has a constant condition
		if (Optimizer<dim>::log_verbose_ && current_iter % interval == 0)
		{
			// Print Iters,value,violation,threshold,feasible_particles to the history
			history_out << current_iter << ","
						<< swarm_[global_best_index_].get_best_value() << ","
						<< swarm_[global_best_index_].get_best_constraint_violation() << ","
						<< violation_threshold_ << ","
						<< feasible_particles <<  std::endl;
			// Print all the particles position and flag the best one to the simulation stream
			for (int i = 0; i < swarm_size_; i++)
			{
				simulation_out << current_iter << ",";
				for (int d = 0; d < dim; d++)
					simulation_out << swarm_[i].get_best_position()[d] << ",";
				simulation_out << ((i == global_best_index_) ? 1 : 0) << std::endl;
			}
		}

		// Update the current iteration
		current_iter++;
	}
}

template <std::size_t dim>
void SASPSO<dim>::print_results(std::ostream &out) const
{
	out << std::setprecision(20) << "Best value found: " << swarm_[global_best_index_].get_best_value() << std::endl;
	out << "Best position found: (";
	for (std::size_t i = 0; i < dim; ++i)
		out << swarm_[global_best_index_].get_best_position()[i] << ", ";
	out << "\b\b)" << std::endl;
	out << "Total constraint violation obtained: " << swarm_[global_best_index_].get_best_constraint_violation() << std::endl;
}

template <std::size_t dim>
bool SASPSO<dim>::is_feasible_solution() const
{
	return swarm_[global_best_index_].get_best_constraint_violation() <= tol_;
}