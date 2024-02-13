#include <memory>

#include "SASPSO/SASPSO.hpp"
#include "TestProblems.hpp"

using namespace type_traits;

int main() {
	auto problem = TestProblems::create_problem<2>(TestProblems::TOWNSEND);

	std::cout << "--- Serial optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<2>> opt = std::make_unique<SASPSO<2>>(problem, 200, 2000);
	opt->initialize();
	opt->optimize();
	opt->print_results();
	std::cout << "Absolute error: " << std::abs(TestProblems::get_exact_value<2>(TestProblems::TOWNSEND) - opt->get_global_best_value()) << std::endl;

	std::cout << "--- Parallel optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<2>> opt_p = std::make_unique<SASPSO<2>>(problem, 200, 2000);
	auto saspso_ptr = dynamic_cast<SASPSO<2>*>(opt_p.get());
    if (!saspso_ptr)
	{
		std::cerr << "Error: dynamic_cast failed" << std::endl;
		return 1;
	}
	saspso_ptr->initialize();
	saspso_ptr->optimize_parallel();
	saspso_ptr->print_results();
	std::cout << "Absolute error: " << std::abs(TestProblems::get_exact_value<2>(TestProblems::TOWNSEND) - opt->get_global_best_value()) << std::endl;

	return 0;
}