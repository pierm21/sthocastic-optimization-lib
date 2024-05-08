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

	// Initialize the bees in the colony
	for (std::size_t i = 0; i < colony_size_; ++i)
	{
		// Create and initialize the bee
		colony_.emplace_back(problem, generator, i);
		colony_[i].initialize();
	}

	global_best_position_ = colony_[0].get_position();
	global_best_value_ = colony_[0].get_value();
	global_best_constraint_violation_ = colony_[0].get_constraint_violation();

	// Find the best bee in the colony
	for (size_t i = 1; i < colony_size_; ++i)
	{
		if (colony_[i].feasibility_rule(global_best_value_, global_best_constraint_violation_))
		{
			global_best_position_ = colony_[i].get_position();
			global_best_value_ = colony_[i].get_value();
			global_best_constraint_violation_ = colony_[i].get_constraint_violation();
		}
	}
}

template <std::size_t dim>
void ABC<dim>::optimize()
{
	int current_iter = 0;

	std::ofstream file_out;
	file_out.open("../output/ABC_output.txt");

	std::uniform_real_distribution<double> distr(0, 1.0);
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	//  Outer optimization loop over all the iterations
	while (current_iter < max_iter_)
	{
		////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////EMPLOYER BEE PHASE///////////////////////////////
		////////////////////////////////////////////////////////////////////////////////

		// Update the position of the bee and prepare data necessary for the onlooker bee phase
		double total_fitness_value = 0;
		double total_constraint_violation = 0;
		for (size_t i = 0; i < colony_size_; ++i)
		{
			// Update the particle
			colony_[i].update_position(MR_, colony_);
			// Prepare the roulette wheel
			total_fitness_value += colony_[i].compute_fitness_value();
			total_constraint_violation += colony_[i].get_constraint_violation();
		}

		////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////ONLOOKER BEE PHASE///////////////////////////////
		////////////////////////////////////////////////////////////////////////////////

		// Compute the probability of each position to be selected by the onlooker bees
		for (size_t i = 0; i < colony_size_; ++i)
		{
			colony_[i].compute_probability(total_fitness_value, total_constraint_violation);
		}

		// Select, according to the fitness probability, the bees that will be updated
		// In other words the position onlooker bees choose to go to and succesively update.
		unsigned int i = 0, t = 0;
		while (t < colony_size_)
		{

			if (distr(generator) < colony_[i].get_probability())
			{
				// file_out << "updated: "<< i << std::endl;
				t++;
				colony_[i].update_position(MR_, colony_);
			}
			i = (i + 1) % colony_size_;
		}

		////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////SCOUT BEE PHASE//////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////

		// If the number of iterations is a multiple of SPP, reset the bees whose position has not
		// been improved since a time that exceeds the limit value
		if (current_iter % SPP_ == 0)
		{
			size_t max_failure_index = 0;
			for (size_t i = 1; i < colony_size_; ++i)
			{
				if (colony_[i].get_failure_counter() >= colony_[max_failure_index].get_failure_counter())
				{
					max_failure_index = i;
				}
			}

			if (colony_[max_failure_index].get_failure_counter() > limit_)
			{
				colony_[max_failure_index].initialize();
			}
		}

		////////////////////////////////////////////////////////////////////////////////

		// Find the best bee in the colony
		for (size_t i = 0; i < colony_size_; ++i)
		{
			if (colony_[i].feasibility_rule(global_best_value_, global_best_constraint_violation_))
			{
				global_best_position_ = colony_[i].get_position();
				global_best_value_ = colony_[i].get_value();
				global_best_constraint_violation_ = colony_[i].get_constraint_violation();
			}
		}

		// Update the current iteration
		current_iter++;
		//file_out << "Best value: " << global_best_value_ << "       Violation: "<< global_best_constraint_violation_<< std::endl;
	}
}

template <std::size_t dim>
void ABC<dim>::print_results(std::ostream &out)
{
	out << std::setprecision(20) << "Best value: " << global_best_value_ << std::endl;
	out << "Best position: (";
	for (std::size_t i = 0; i < dim; ++i)
		out << global_best_position_[i] << ", ";
	out << "\b\b)" << std::endl;
	out << "Total constraint violation: " << global_best_constraint_violation_ << std::endl;
}

template <std::size_t dim>
void ABC<dim>::initialize_parallel() {}

template <std::size_t dim>
void ABC<dim>::optimize_parallel() {}

template <std::size_t dim>
void ABC<dim>::print_initizalization(std::ostream &out)
{
	int b = 0;
	for (int j = 0; j < colony_size_; j++)
	{
		std::cout << "Bee " << b++ << ": " << std::endl;
		for (int i = 0; i < dim; i++)
		{
			std::cout << colony_[j].get_position()[i] << " " << std::endl;
		}
		std::cout << std::endl;
	}
}
