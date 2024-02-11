#pragma once

#include <iostream>

#include "TypeTraits.hpp"
#include "Problem.hpp"

template <size_t dim>
class Optimizer {
protected:
	Problem<dim> problem_;

public:
	Optimizer() = delete;
	Optimizer(const Problem<dim>& problem)
		: problem_(problem) {}
	virtual ~Optimizer() = default;
	virtual void initialize() = 0;
	virtual void optimize() = 0;
	virtual void print_results(std::ostream& out = std::cout) = 0;

	// getters
	virtual double get_global_best_value() = 0;
	virtual const RealVector<dim>& get_global_best_position() = 0;
	virtual bool is_feasible_solution() = 0;
};