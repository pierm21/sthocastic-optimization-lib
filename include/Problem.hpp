#pragma once

#include "TypeTraits.hpp"

template <size_t dim>
class Problem
{
private:
	RealFunction<dim> fitness_function_;
	RealVector<dim> lower_bounds_;
	RealVector<dim> upper_bounds_;
	std::vector<RealFunction<dim>> equality_constraints_;
	std::vector<RealFunction<dim>> inequality_constraints_;

public:
	Problem(const RealFunction<dim> &fitness_function, const RealVector<dim> &lower_bounds, const RealVector<dim> &upper_bounds)
		: fitness_function_(fitness_function),
		  lower_bounds_(lower_bounds),
		  upper_bounds_(upper_bounds){};

	~Problem() = default;

	const RealFunction<dim> &get_fitness_function() const { return fitness_function_; }
	const RealVector<dim> &get_lower_bounds() const { return lower_bounds_; }
	const RealVector<dim> &get_upper_bounds() const { return upper_bounds_; }
	double get_lower_bound(size_t i) const { return lower_bounds_[i]; }
	double get_upper_bound(size_t i) const { return upper_bounds_[i]; }
	const std::vector<RealFunction<dim>> &get_equality_constraints() const { return equality_constraints_; }
	const std::vector<RealFunction<dim>> &get_inequality_constraints() const { return inequality_constraints_; }

	void add_equality_constraint(RealFunction<dim> &constraint) { equality_constraints_.push_back(constraint); }
	void add_inequality_constraint(RealFunction<dim> &constraint) { inequality_constraints_.push_back(constraint); }

	 bool has_constraints() const { return !equality_constraints_.empty() || !inequality_constraints_.empty(); }
};