#pragma once

#include <iostream>

#include "TypeTraits.hpp"
#include "Problem.hpp"

template <size_t dim>
class Optimizer {
protected:
	Problem<dim> problem_;
	double global_best_value_;
	RealVector<dim> global_best_position_;

public:
	Optimizer() = delete;
	Optimizer(const Problem<dim>& problem)
		: problem_(problem) {}
	virtual void initialize() = 0;
	virtual void optimize() = 0;
	virtual void print_results(std::ostream& out = std::cout) = 0;
};