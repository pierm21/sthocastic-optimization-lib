#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <filesystem>

#include "ABC/ABC.hpp"
#include "TestProblems.hpp"
#include "OptimizerFactory.hpp"

using namespace type_traits;
namespace fs = std::filesystem;

/*#define dimension 2
#define test_problem TestProblems::GOMEZ_LEVY
#define problem_name "GOMEZ_LEVY"*/

/*#define dimension 2
#define test_problem TestProblems::TOWNSEND
#define problem_name "TOWNSED"*/

/*#define dimension 8
#define test_problem TestProblems::G10
#define problem_name "G10"*/

/*#define dimension 10
#define test_problem TestProblems::G7
#define problem_name "G7"*/


int optimize()
{
	constexpr int log_interval = 50;
	std::vector<double> best_values;
	int iter = 6000;
	int particles = 20;
	double tol = 1e-16;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	// Preliminary informations to std out
	std::cout << "Parallel adaptive optimization" << std::endl;
	std::cout << "Logs in /output/saspso_optimize.csv" << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Initialize history
	std::vector<double> history, violation, feasible, threshold;		

	for (int i = 0; i < 30; i++)
	{
	std::unique_ptr<Optimizer<dimension>> opt_p = std::make_unique<ABC<dimension>>(problem, particles, iter);//, 1e-10);
	auto abc_ptr = dynamic_cast<ABC<dimension> *>(opt_p.get());
	if (!abc_ptr)
	{
		std::cerr << "Error: dynamic_cast failed" << std::endl;
		return 1;
	}

	std::ofstream file_out;
	file_out.open("../output/abc_optimize.csv");
	if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	// Write comments and header
	file_out << "# Fitness, constraints violation and feasible bees over iterations" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Bees: " << particles << std::endl;
	file_out << "# Tolerance: " << tol << std::endl;
	file_out << "# Problem: " << problem_name << std::endl;

	file_out << "Iters,value,violation,threshold,feasible_bees" << std::endl;

	//saspso_ptr->initialize_parallel();
	//saspso_ptr->optimize_parallel(history, violation, feasible, threshold, log_interval);
    abc_ptr->initialize();
    //abc_ptr->print_initizalization();
    abc_ptr->optimize();
    //abc_ptr->optimize(/*history, violation, feasible, threshold, log_interval*/);

	// Get the exact global minimum
	double exact_value = TestProblems::get_exact_value<dimension>(test_problem);
    std::cout << "Exact value: " << exact_value << std::endl;

    // Get the minimum value found
    double found_value = abc_ptr->get_global_best_value();
    std::cout << "Found value: " << found_value << std::endl;
	std::cout << std::endl;
	best_values.push_back(found_value);
	}

	// Get the minimum value found
	double found_value = best_values[0];
	for (int i = 1; i < best_values.size(); i++)
	{
		if (best_values[i] < found_value)
			found_value = best_values[i];
	}
	std::cout << "Minimum Found value: " << found_value << std::endl;

	// Compute mean and standard deviation
	double mean = 0;
	for (int i = 0; i < best_values.size(); i++)
	{
		mean += best_values[i];
	}
	mean /= best_values.size();
	std::cout << "Mean: " << mean << std::endl;

	double std_dev = 0;
	for (int i = 0; i < best_values.size(); i++)
	{
		std_dev += (best_values[i] - mean) * (best_values[i] - mean);
	}
	std_dev = std::sqrt(std_dev / best_values.size());
	std::cout << "Standard deviation: " << std_dev << std::endl;

	// Print the final error
	//std::cout << std::endl << "Absolute error: " << std::abs(history.back() - exact_value) << std::endl;
	//std::cout << "Relative error: " << std::abs(history.back() - exact_value) / exact_value << std::endl;
    //std::cout << std::endl << "Absolute error: " << std::abs(found_value - exact_value) << std::endl;
	//std::cout << "Relative error: " << std::abs(found_value - exact_value) / exact_value << std::endl;

/*	// Store on file
	for (int i = 0; i < history.size(); i++)
	{
		file_out << i * log_interval << ",";
		file_out << history[i] << ",";
		file_out << violation[i] << ",";
		file_out << threshold[i] << ",";
		file_out << feasible[i] << std::endl;
	}
	file_out.close();
*/
	return 0;
}

int main(int argc, char **argv)
{
	// Check the number of arguments
	if (argc != 2)
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for SASPSO algorithm: static_adaptive, serial_parallel, time_numparticles" << std::endl;
		return -1;
	}
	// Create if it not exist the output directory
	fs::create_directory("../output");

	// Get from command line the required test
	std::string test = argv[1];
	/*if (test == "static_adaptive")
		static_adaptive_test();
	else if (test == "serial_parallel")
		serial_parallel_test();
	else if (test == "time_numparticles")
		time_numparticles_test();
	else*/
    if (test == "optimize")
		optimize();
	else
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for SASPSO algorithm: static_adaptive, serial_parallel, time_numparticles" << std::endl;
		return -1;
	}

	return 0;
}