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

#define dimension 8
#define test_problem TestProblems::G10
#define problem_name "G10"

int serial_parallel_test()
{
	// Initialize problem and solver parameters
	int iter = 5000;
	int particles = 300;
	double tol = 1e-16;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Test the serial version
	std::cout << std::endl;
	std::cout << "--- Serial optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<dimension>> opt = std::make_unique<SASPSO<dimension>>(problem, particles, iter, tol);
	opt->initialize();
	auto t1 = std::chrono::high_resolution_clock::now();
	opt->optimize();
	auto t2 = std::chrono::high_resolution_clock::now();
	opt->print_results();

	std::cout << std::setprecision(20) << "Absolute error: " << std::abs(TestProblems::get_exact_value<dimension>(test_problem) - opt->get_global_best_value()) << std::endl;
	std::cout << std::setprecision(20) << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - opt->get_global_best_position()).norm() << std::endl;
	std::cout << std::setprecision(20) << "Elapsed optimization time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;

	// Test the parallel version
	std::cout << std::endl;
	std::cout << "--- Parallel optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<dimension>> opt_p = std::make_unique<SASPSO<dimension>>(problem, particles, iter, 1e-10);
	auto saspso_ptr = dynamic_cast<SASPSO<dimension> *>(opt_p.get());
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
	std::cout << std::setprecision(20) << "Absolute error: " << std::abs(TestProblems::get_exact_value<dimension>(test_problem) - opt->get_global_best_value()) << std::endl;
	std::cout << std::setprecision(20) << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - opt->get_global_best_position()).norm() << std::endl;
	std::cout << std::setprecision(20) << "Elapsed optimization time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;
	return 0;
}

int static_adaptive_test()
{
	constexpr int log_interval = 10;

	int iter = 2300;
	int particles = 200;
	double tol = 1e-16;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	// Preliminary informations to std out
	std::cout << "Error as function of performed iteration test (for static and adaptive)." << std::endl;
	std::cout << "Logs in /output/saspso_static_adaptive.csv every " << log_interval << " iterations." << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Initialize the file
	std::ofstream file_out;
	file_out.open("../output/saspso_static_adaptive.csv");
	if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	// Write comments and header
	file_out << "# Error and constraints violation over iterations (static vs adaptive)" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Particles: " << particles << std::endl;
	file_out << "# Tolerance: " << tol << std::endl;
	file_out << "# Problem: " << problem_name << std::endl;

	file_out << "Iters,static_err,static_viol,static_p_err,static_p_viol,adaptive_err,adaptive_viol,adaptive_p_err,adaptive_p_viol" << std::endl;

	// Initialize history
	std::vector<double> history_s, violation_s;
	std::vector<double> history_d, violation_d;
	std::vector<double> history_s_p, violation_s_p;
	std::vector<double> history_d_p, violation_d_p;
	std::vector<double> dummy;

	//TODO: are pointers necessary? Couldn't we just use the object directly?
	// Optimize the problems. We need to use a pointer for the specialized class
	std::unique_ptr<SASPSO<dimension>> opt = std::make_unique<SASPSO<dimension>>(problem, particles, iter, tol);
	opt->initialize();
	opt->optimize(history_d, violation_d, log_interval);

	opt = std::make_unique<SASPSO<dimension>>(problem, particles, iter, tol);
	opt->initialize_parallel();
	opt->optimize_parallel(history_d_p, violation_d_p, dummy, dummy, log_interval);

	std::unique_ptr<SASPSO<dimension>> opt1 = std::make_unique<SASPSO<dimension>>(problem, particles, iter, tol, 0.9, 0.9, 1.2, 1.2, 0.3, 0.3);
	opt1->initialize();
	opt1->optimize(history_s, violation_s, log_interval);

	opt1 = std::make_unique<SASPSO<dimension>>(problem, particles, iter, tol, 0.9, 0.9, 1.2, 1.2, 0.3, 0.3);
	opt1->initialize_parallel();
	opt1->optimize_parallel(history_s_p, violation_s_p, dummy, dummy, log_interval);

	// Print the history vectors size
	std::cout << "Adaptive history size: " << history_d.size() << std::endl;
	std::cout << "Adaptive parallel history size: " << history_d_p.size() << std::endl;
	std::cout << "Static history size: " << history_s.size() << std::endl;
	std::cout << "Static parallel history size: " << history_s_p.size() << std::endl;

	// Get the exact global minimum
	double exact_value = TestProblems::get_exact_value<dimension>(test_problem);

	// Write to file the error values
	for (int i = 0; i < history_s.size(); i++)
	{
		file_out << i * log_interval << ",";
		file_out << std::abs(history_s[i] - exact_value) << ",";
		file_out << violation_s[i] << ",";
		file_out << std::abs(history_s_p[i] - exact_value) << ",";
		file_out << violation_s_p[i] << ",";
		file_out << std::abs(history_d[i] - exact_value) << ",";
		file_out << violation_d[i] << ",";
		file_out << std::abs(history_d_p[i] - exact_value) << ",";
		file_out << violation_d_p[i] << std::endl;
	}

	// Close the file
	file_out.close();
	return 0;
}

// test the time as function of the number of particles
int time_numparticles_test()
{
	int iter = 2000;
	double tol = 1e-10;
	int max_particles = 2000;
	constexpr int log_interval = 50;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	// Preliminary informations to std out
	std::cout << "Time and Speedup as function of swarm size." << std::endl;
	std::cout << "Logs in /output/saspso_time_numparticles.csv" << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Log interval: " << log_interval << std::endl;

	// Get the number of omp threads
	int tot_threads = 0;
#pragma omp parallel
	{
#pragma omp single
		tot_threads = omp_get_num_threads();
	}

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
	file_out << "# Problem: " << problem_name << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Threads: " << tot_threads << std::endl;
	file_out << "Num_particles,Serial_time,Parallel_time,Speedup" << std::endl;

	std::cout << "Starting test from 1 to " << max_particles << " particles" << std::endl;
	std::cout << "Logging every " << log_interval << " iterations" << std::endl;

	for (int i = log_interval; i <= max_particles; i += log_interval)
	{
		// Print progress to stdout
		if ((i) % (log_interval) == 0)
			std::cout << "Starting test with " << i << " particle(s)" << std::endl;
		// Optimize the problem serially
		std::unique_ptr<SASPSO<dimension>> opt = std::make_unique<SASPSO<dimension>>(problem, i, iter, tol);
		opt->initialize();
		auto t1 = std::chrono::high_resolution_clock::now();
		opt->optimize();
		auto t2 = std::chrono::high_resolution_clock::now();
		// Optimize parallel
		opt = std::make_unique<SASPSO<dimension>>(problem, i, iter, tol);
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

int optimize()
{
	constexpr int log_interval = 50;

	int iter = 20000;
	int particles = 5000;
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

	std::unique_ptr<Optimizer<dimension>> opt_p = std::make_unique<SASPSO<dimension>>(problem, particles, iter, 1e-10);
	auto saspso_ptr = dynamic_cast<SASPSO<dimension> *>(opt_p.get());
	if (!saspso_ptr)
	{
		std::cerr << "Error: dynamic_cast failed" << std::endl;
		return 1;
	}

	std::ofstream file_out;
	file_out.open("../output/saspso_optimize.csv");
	if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	// Write comments and header
	file_out << "# Fitness, constraints violation and feasible particles over iterations" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Particles: " << particles << std::endl;
	file_out << "# Tolerance: " << tol << std::endl;
	file_out << "# Problem: " << problem_name << std::endl;

	file_out << "Iters,value,violation,threshold,feasible_particles" << std::endl;

	saspso_ptr->initialize_parallel();
	saspso_ptr->optimize_parallel(history, violation, feasible, threshold, log_interval);

	// Get the exact global minimum
	double exact_value = TestProblems::get_exact_value<dimension>(test_problem);

	// Print the final error
	std::cout << std::endl << "Absolute error: " << std::abs(history.back() - exact_value) << std::endl;
	std::cout << "Relative error: " << std::abs(history.back() - exact_value) / exact_value << std::endl;

	// Store on file
	for (int i = 0; i < history.size(); i++)
	{
		file_out << i * log_interval << ",";
		file_out << history[i] << ",";
		file_out << violation[i] << ",";
		file_out << threshold[i] << ",";
		file_out << feasible[i] << std::endl;
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
		std::cout << "Available tests for SASPSO algorithm: static_adaptive, serial_parallel, time_numparticles" << std::endl;
		return -1;
	}
	// Create if it not exist the output directory
	fs::create_directory("../output");

	// Get from command line the required test
	std::string test = argv[1];
	if (test == "static_adaptive")
		static_adaptive_test();
	else if (test == "serial_parallel")
		serial_parallel_test();
	else if (test == "time_numparticles")
		time_numparticles_test();
	else if (test == "optimize")
		optimize();
	else
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for SASPSO algorithm: static_adaptive, serial_parallel, time_numparticles" << std::endl;
		return -1;
	}

	return 0;
}