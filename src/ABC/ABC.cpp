#pragma once

#include "ABC/ABC.hpp"
#include <omp.h>
#include <iomanip>
#include "TestProblems.hpp"
#include "OptimizerFactory.hpp"

template <std::size_t dim>
void ABC<dim>::initialize()
{
	// Instantiate the marsenne twister
	std::random_device rand_dev;
	std::shared_ptr<std::mt19937> generator = std::make_shared<std::mt19937>(rand_dev());
	// Create a shared pointer to the problem
	std::shared_ptr<Problem<dim>> problem = std::make_shared<Problem<dim>>(Optimizer<dim>::problem_);

    // Initialize the array of total constraint violations
	std::vector<double> total_violations;

    // Initialize the bees in the colony
	for (std::size_t i = 1; i < colony_size_; ++i)
	{
		// Create and initialize the bee
		colony_.emplace_back(problem, generator);
		colony_[i].initialize();
		// Add the constraint violation to the array
		total_violations.push_back(colony_[i].get_constraint_violation());
	}

    // Initialize the violation threshold as the median of the total violations
	std::sort(total_violations.begin(), total_violations.end());
	if (colony_size_ % 2 == 0)
		violation_threshold_ = (total_violations[colony_size_ / 2] + total_violations[colony_size_ / 2 - 1]) / 2;
	else
		violation_threshold_ = total_violations[colony_size_ / 2];
}


template <std::size_t dim>
void ABC<dim>::optimize()
{
	int current_iter = 0;

	// Outer optimization loop over all the iterations
	while (current_iter < max_iter_)
	{
		// Reset the number of feasible bees for the current iteration
		int feasible_bees_ = 0;

        ////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////EMPLOYER BEE PHASE///////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

		for (size_t i = 0; i < colony_size_; ++i)
		{
			// Update the particle
			colony_[i].update_position(MR_, violation_threshold_, &colony_, tol_);
			// Check if the particle is feasible
			if (colony_[i].get_constraint_violation() <= violation_threshold_)
				feasible_bees_++;
		}

		// Update the violation threshold according to the proportion of feasible particles
		violation_threshold_ = violation_threshold_ * (1.0 - (feasible_bees_ / (double)colony_size_));
		violation_threshold_ = violation_threshold_ < tol_ ? 0 : violation_threshold_;

		
        ////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////ONLOOKER BEE PHASE///////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        
        


        // Update the current iteration
		current_iter++; 
	}
}
