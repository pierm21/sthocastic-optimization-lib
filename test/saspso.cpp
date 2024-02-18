#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <filesystem>

#include "SASPSO/SASPSO.hpp"
#include "TestProblems.hpp"
#include "OptimizerFactory.hpp"

using namespace type_traits;
namespace fs = std::filesystem;

#define dimension 2

int serial_parallel_test()
{
	// Initialize problem and solver parameters
	int iter = 5000;
	int particles = 300;
	auto problem = TestProblems::create_problem<2>(TestProblems::GOMEZ_LEVY);

	std::cout << "Problem: Gomez Levy" << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Test the serial version
	std::cout << "--- Serial optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<2>> opt = std::make_unique<SASPSO<2>>(problem, particles, iter, 1e-10);
	opt->initialize();
	auto t1 = std::chrono::high_resolution_clock::now();
	opt->optimize();
	auto t2 = std::chrono::high_resolution_clock::now();
	opt->print_results();
	std::cout << "Absolute error: " << std::abs(TestProblems::get_exact_value<2>(TestProblems::GOMEZ_LEVY) - opt->get_global_best_value()) << std::endl;
	std::cout << "Absolute distance: " << (TestProblems::get_exact_position<2>(TestProblems::GOMEZ_LEVY) - opt->get_global_best_position()).norm() << std::endl;
	std::cout << "Elapsed optimization time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;

	// Yest the parallel version
	std::cout << "--- Parallel optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<2>> opt_p = std::make_unique<SASPSO<2>>(problem, particles, iter, 1e-10);
	auto saspso_ptr = dynamic_cast<SASPSO<2> *>(opt_p.get());
	if (!saspso_ptr)
	{
		std::cerr << "Error: dynamic_cast failed" << std::endl;
		return 1;
	}
	saspso_ptr->initialize_parallel();
	t1 = std::chrono::high_resolution_clock::now();
	saspso_ptr->optimize_parallel();
	t2 = std::chrono::high_resolution_clock::now();
	saspso_ptr->print_results();
	std::cout << "Absolute error: " << std::abs(TestProblems::get_exact_value<2>(TestProblems::GOMEZ_LEVY) - opt->get_global_best_value()) << std::endl;
	std::cout << "Absolute distance: " << (TestProblems::get_exact_position<2>(TestProblems::GOMEZ_LEVY) - opt->get_global_best_position()).norm() << std::endl;
	std::cout << "Elapsed optimization time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;
	return 0;
}

int error_iterations_test()
{
	constexpr int log_interval = 10;

	int iter = 10000;
	int particles = 300;
	double tol = 1e-16;
	auto gomez_levy = TestProblems::create_problem<dimension>(TestProblems::GOMEZ_LEVY);
	auto townsend = TestProblems::create_problem<dimension>(TestProblems::TOWNSEND);

	// Preliminary informations to std out
	std::cout << "Error as function of performed iteration test." << std::endl;
	std::cout << "Logs in /output/saspso_error_iterations.csv every " << log_interval << " iterations." << std::endl;

	std::cout << "Problem: Townsend, Gomez-Levy" << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Initialize the file
	std::ofstream file_out;
	file_out.open("../output/saspso_error_iterations.csv");
	if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	// Write comments and header
	file_out << "# Error as function of performed iterations for Townsend and Gomez-Levy problems" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Particles: " << particles << std::endl;
	file_out << "# Tolerance: " << tol << std::endl;

	file_out << "Iterations,Townsend,Townsend_parallel,Gomez-Levy,Gomez-Levy_parallel" << std::endl;

	// Initialize history and log_interval
	std::vector<double> history_townsend;
	std::vector<double> history_gomez_levy;
	std::vector<double> history_townsend_p;
	std::vector<double> history_gomez_levy_p;

	// Optimize the problems. We need to use a pointer for the specialized class
	std::unique_ptr<SASPSO<2>> opt = std::make_unique<SASPSO<2>>(gomez_levy, particles, iter, tol);
	opt->initialize();
	opt->optimize(history_gomez_levy, log_interval);

	opt = std::make_unique<SASPSO<2>>(gomez_levy, particles, iter, tol);
	opt->initialize_parallel();
	opt->optimize_parallel(history_gomez_levy_p, log_interval);

	opt = std::make_unique<SASPSO<2>>(townsend, particles, iter, tol);
	opt->initialize();
	opt->optimize(history_townsend, log_interval);

	opt = std::make_unique<SASPSO<2>>(townsend, particles, iter, tol);
	opt->initialize_parallel();
	opt->optimize_parallel(history_townsend_p, log_interval);

	// Print the history vectors size
	std::cout << "Townsend history size: " << history_townsend.size() << std::endl;
	std::cout << "Townsend parallel history size: " << history_townsend_p.size() << std::endl;
	std::cout << "Gomez-Levy history size: " << history_gomez_levy.size() << std::endl;
	std::cout << "Gomez-Levy parallel history size: " << history_gomez_levy_p.size() << std::endl;

	// Get the exact global minimum
	double exact_townsend = TestProblems::get_exact_value<dimension>(TestProblems::TOWNSEND);
	double exact_gomez_levy = TestProblems::get_exact_value<dimension>(TestProblems::GOMEZ_LEVY);

	// Write to file the error values
	for (int i = 0; i < history_townsend.size(); i++)
	{
		file_out << i * log_interval << ",";
		file_out << std::abs(history_townsend[i] - exact_townsend) << ",";
		file_out << std::abs(history_townsend_p[i] - exact_townsend) << ",";
		file_out << std::abs(history_gomez_levy[i] - exact_gomez_levy) << ",";
		file_out << std::abs(history_gomez_levy_p[i] - exact_gomez_levy) << std::endl;
	}

	// Close the file
	file_out.close();
	return 0;
}

// test the time as function of the number of particles
int time_numparticles_test()
{
	int iter = 3000;
	double tol = 1e-10;
	int max_particles = 500;
	constexpr int log_interval = 10;
	auto gomez_levy = TestProblems::create_problem<dimension>(TestProblems::GOMEZ_LEVY);

	// Preliminary informations to std out
	std::cout << "Time and Speedup as function of swarm size." << std::endl;
	std::cout << "Logs in /output/saspso_time_numparticles.csv" << std::endl;

	std::cout << "Problem: Townsend" << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Log interval: " << log_interval << std::endl;

	// Initialize the file
	std::ofstream file_out;
	file_out.open("../output/saspso_time_numparticles.csv");
	if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}
	// Write comments and header to file
	file_out << "# Execution time as function of the swarm size" << std::endl;
	file_out << "# Problem: Townsend" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "Num_particles,Serial_time,Parallel_time,Speedup" << std::endl;

	std::cout << "Starting test from 1 to " << max_particles << " particles" << std::endl;
	std::cout << "Logging every " << log_interval << " iterations" << std::endl;

	for (int i = 1; i <= max_particles; i += log_interval)
	{
		// Print progress to stdout
		if ((i - 1) % (log_interval) == 0)
			std::cout << "Starting test with " << i << " particle(s)" << std::endl;
		// Optimize the problem serially
		std::unique_ptr<SASPSO<2>> opt = std::make_unique<SASPSO<2>>(gomez_levy, i, iter, tol);
		opt->initialize();
		auto t1 = std::chrono::high_resolution_clock::now();
		opt->optimize();
		auto t2 = std::chrono::high_resolution_clock::now();
		// Optimize in parallel
		opt = std::make_unique<SASPSO<2>>(gomez_levy, i, iter, tol);
		opt->initialize_parallel();
		auto t3 = std::chrono::high_resolution_clock::now();
		opt->optimize_parallel();
		auto t4 = std::chrono::high_resolution_clock::now();
		// Write data to file
		file_out << i << ",";
		file_out << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << ",";
		file_out << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << ",";
		file_out << double(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count())
			/ std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << std::endl;
	}
	file_out.close();
	return 0;
}

int main(int argc, char **argv)
{
	// Check the number of arguments
	if (argc != 2)
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for SASPSO algorithm: error_iterations, serial_parallel, time_numparticles" << std::endl;
		return -1;
	}
	// Create if it not exist the output directory
	fs::create_directory("../output");

	// Get from command line the required test
	std::string test = argv[1];
	if (test == "error_iterations")
		error_iterations_test();
	else if (test == "serial_parallel")
		serial_parallel_test();
	else if (test == "time_numparticles")
		time_numparticles_test();
	else
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for SASPSO algorithm: error_iterations, serial_parallel, time_numparticles" << std::endl;
		return -1;
	}
	return 0;
}