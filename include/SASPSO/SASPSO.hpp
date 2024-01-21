#pragma once

#include "Optimizer.hpp"

template <std::size_t dim>
class SASPSO : public Optimizer<dim> {

private:
	int swarm_size_;
	double c1_, c2_, w_;

public:
	SASPSO(const Problem<dim>& problem, int max_iter, int swarm_size, double tol, double c1, double c2, double w)
		: Optimizer<dim>(problem), swarm_size_(swarm_size), c1_(c1), c2_(c2), w_(w) {}

	void initialize() override;
	void optimize() override;
	void print_results(std::ostream& out = std::cout) override;

};

// from here starts the SPSPSO.cpp file

template <std::size_t dim>
void SASPSO<dim>::initialize() {
	std::cout << "SASPSO::initialize()" << w_ << std::endl;
	SASPSO::problem_.has_constraints() ? std::cout << "Constraints are present" << std::endl : std::cout << "No constraints" << std::endl;
}

template <std::size_t dim>
inline void SASPSO<dim>::optimize()
{
	std::cout << "SASPSO::optimize()" << std::endl;
}

template <std::size_t dim>
inline void SASPSO<dim>::print_results(std::ostream& out)
{
	std::cout << "SASPSO::print_results()" << std::endl;
}
