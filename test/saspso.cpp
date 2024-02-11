#include <memory>

#include "SASPSO/SASPSO.hpp"
#include "TestProblems.hpp"

using namespace type_traits;

int main() {
	auto problem = TestProblems::create_problem<2>(TestProblems::TOWNSEND);

	// particle testing
	auto problem_ptr = std::make_shared<Problem<2>>(problem);
	std::random_device rand_dev;
	auto gen = std::make_shared<std::mt19937>(rand_dev());
	Particle<2> particle(problem_ptr, gen, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0);
	particle.initialize();
	particle.print();
	particle.update(RealVector<2>({ 0.0, 0.0 }), 1, 100);
	particle.print();

	// optimizer testing
	std::cout << "--- Optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<2>> opt = std::make_unique<SASPSO<2>>(problem, 300, 5000, 1e-8);
	opt->initialize();
	opt->optimize();
	opt->print_results();
	std::cout << "Absolute error: " << std::abs(TestProblems::get_exact_value<2>(TestProblems::TOWNSEND) - opt->get_global_best_value()) << std::endl;

	return 0;
}