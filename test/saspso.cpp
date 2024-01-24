#include <memory>

#include "SASPSO/SASPSO.hpp"
#include "TestProblems.hpp"

using namespace type_traits;

// test function
// f(x) = x1^2 + x2^2
double f(const RealVector<2>& x) {
	return x[0] * x[0] + x[1] * x[1];
}

int main() {
	auto problem = TestProblems::create_problem<2>(TestProblems::TOWNSEND);
	auto problem_ptr = std::make_shared<Problem<2>>(problem);

	// particle testing
	std::random_device rand_dev;
	auto gen = std::make_shared<std::mt19937>(rand_dev());
	Particle<2> particle(problem_ptr, gen, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0);
	particle.initialize();
	particle.print();
	particle.update(RealVector<2>({ 0.0, 0.0 }), 1, 100);
	particle.print();

	// optimizer testing
	std::unique_ptr<Optimizer<2>> opt = std::make_unique<SASPSO<2>>(problem, 100, 2000, 1e-6);
	opt->initialize();
	opt->optimize();
}