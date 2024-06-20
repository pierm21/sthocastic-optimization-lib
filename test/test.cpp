#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <filesystem>
#include <mpi.h>

#include "TestProblems.hpp"
#include "OptimizerFactory.hpp"

using namespace type_traits;
namespace fs = std::filesystem;

/*#define dimension 2
#define test_problem TestProblems::GOMEZ_LEVY
#define problem_name "GOMEZ_LEVY"
*/

/*#define dimension 2
#define test_problem TestProblems::TOWNSEND
#define problem_name "TOWNSEND"*/

/*#define dimension 8
#define test_problem TestProblems::G10
#define problem_name "G10"*/

#define dimension 30
#define test_problem TestProblems::GRIEWANK
#define problem_name "GRIEWANK"

/*#define dimension 10
#define test_problem TestProblems::G7
#define problem_name "G7"*/

/**
 * @brief Optimize the given problem using the specified algorithm and log the results to the output files
 * Logging is made both for history and simulation data
 * @param algorithm_name the enum of the algorithm to be used for optimization
 * @return int 0 if the optimization is successful, -1 otherwise
 */
int optimize(const typename OptimizerFactory<dimension>::OptimizerName &algorithm_name)
{
	constexpr int log_interval = 50;

	int iter = 6000;
	int particles = 20;
	auto problem = TestProblems::create_problem<dimension>(test_problem);
	std::string alg_name_str = OptimizerFactory<dimension>::get_string_name(algorithm_name);

	// Preliminary informations to std out
	std::cout << "Serial optimization for " << alg_name_str << std::endl;
	std::cout << "Logs in /output/" + alg_name_str + "_optimize.csv" << std::endl;
	std::cout << "Simulation data in /output/" + alg_name_str + "_simulation.csv" << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	std::unique_ptr<Optimizer<dimension>> opt = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter);
	opt->set_log_verbose(true);

	std::ofstream history_out;
	history_out.open("../output/optimize_" + alg_name_str + "_" + std::to_string(1) + ".csv");
	if (!history_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}
	std::ofstream simulation_out;
	simulation_out.open("../output/simulation_" + alg_name_str + "_" + std::to_string(1) + ".csv");
	if (!simulation_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	// Write comments and header
	history_out << "# Fitness, constraints violation and feasible particles over iterations" << std::endl;
	history_out << "# Dimension: " << dimension << std::endl;
	history_out << "# Particles: " << particles << std::endl;
	history_out << "# Problem: " << problem_name << std::endl;
	history_out << "iters,value,violation,feasible_particles" << std::endl;

	simulation_out << "# Data of all the particles position at every interval iterations, the best one is flagged" << std::endl;
	simulation_out << "# Dimension: " << dimension << std::endl;
	simulation_out << "# Particles: " << particles << std::endl;
	simulation_out << "# Problem: " << problem_name << std::endl;
	simulation_out << "iter,";
	for (int i = 0; i < dimension; i++)
		simulation_out << "x" << i << ",";
	simulation_out << "isbest" << std::endl;

	// Optimize the problem storing data to files
	opt->initialize();
	opt->optimize(history_out, simulation_out, log_interval);

	// Get the exact global minimum
	double exact_value = TestProblems::get_exact_value<dimension>(test_problem);

	// Print the final data
	opt->print_results();
	std::cout << "Exact value: " << TestProblems::get_exact_value<dimension>(test_problem) << std::endl;
	std::cout << "Absolute error: " << std::abs(opt->get_global_best_value() - exact_value) << std::endl;
	std::cout << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - opt->get_global_best_position()).norm() << std::endl;

	// Close file streams
	history_out.close();
	simulation_out.close();

	return 0;
}

// test the time as function of the number of particles
int time_numparticles_test(const typename OptimizerFactory<dimension>::OptimizerName &algorithm_name)
{
	int iter = 10000;
	int max_particles = 2048;
	constexpr int log_multiplier = 2;
	constexpr int init_particles = 32;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	std::string alg_name_str = OptimizerFactory<dimension>::get_string_name(algorithm_name);

	// Get the number of omp threads
	int tot_threads = 0;
#pragma omp parallel
	{
#pragma omp single
		tot_threads = omp_get_num_threads();
	}

	// Preliminary informations to std out
	std::cout << "Time and Speedup as function of colony size for " << alg_name_str << std::endl;
	std::cout << "Logs in /output/time_numparticles_" + alg_name_str + "_" + std::to_string(tot_threads) + ".csv" << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Log multiplier: " << log_multiplier << std::endl;


	// Initialize the file
	std::ofstream file_out;
	file_out.open("../output/time_numparticles_" + alg_name_str + "_" + std::to_string(tot_threads) + ".csv");
	if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}
	// Write comments and header to file
	file_out << "# Execution time as function of the colony size" << std::endl;
	file_out << "# Problem: " << problem_name << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Threads: " << tot_threads << std::endl;
	file_out << "Num_particles,Serial_time,Parallel_time,Speedup" << std::endl;

	std::cout << "Starting test from 1 to " << max_particles << " particles" << std::endl;
	std::cout << "Logging every x" << log_multiplier << " iterations" << std::endl;

	for (int i = init_particles; i <= max_particles; i *= log_multiplier)
	{
		// Print progress to stdout
		std::cout << "Starting test with " << i << " particle(s)" << std::endl;
		// Optimize the problem serially
		std::unique_ptr<Optimizer<dimension>> opt = OptimizerFactory<dimension>::create(algorithm_name, problem, i, iter);
		opt->set_log_verbose(false);

		opt->initialize();
		auto t1 = std::chrono::high_resolution_clock::now();
		opt->optimize();
		auto t2 = std::chrono::high_resolution_clock::now();

		// Optimize parallel
		opt = OptimizerFactory<dimension>::create(algorithm_name, problem, i, iter);
		opt->set_log_verbose(false);

		opt->initialize_parallel();
		auto t3 = std::chrono::high_resolution_clock::now();
		opt->optimize_parallel();
		auto t4 = std::chrono::high_resolution_clock::now();

		// Write data to file
		file_out << i << ",";
		file_out << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << ",";
		file_out << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << ",";
		file_out << double(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()) / std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count() << std::endl;
	}
	file_out.close();
	return 0;

}

int static_adaptive_test(const typename OptimizerFactory<dimension>::OptimizerName &algorithm_name)
{ /*
	constexpr int log_interval = 10;

	int iter = 2300;
	int particles = 200;
	double tol = 1e-16;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	std::string alg_name_str = OptimizerFactory<dimension>::get_string_name(algorithm_name);

	// Preliminary informations to std out
	std::cout << "Error as function of performed iteration test (for static and adaptive, serial parallel)." << std::endl;
	std::cout << "Logs in /output/" << alg_name_str << "_static_adaptive.csv every " << log_interval << " iterations." << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Initialize the file
	std::ofstream file_out, temp_1;
	file_out.open("../output/" + alg_name_str + "_static_adaptive.csv");
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

	file_out << "iters,static_err,static_viol,static_p_err,static_p_viol,adaptive_err,adaptive_viol,adaptive_p_err,adaptive_p_viol" << std::endl;

	// Initialize history
	std::vector<double> history_s, violation_s, feasible_s;
	std::vector<double> history_d, violation_d, feasible_d;
	std::vector<double> history_s_p, violation_s_p, feasible_s_p;
	std::vector<double> history_d_p, violation_d_p, feasible_d_p;
	std::vector<double> dummy;

	// Optimize the problems. We need to use a pointer for the specialized class
	std::unique_ptr<Optimizer<dimension>> opt = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter, tol);
	opt->initialize();
	opt->optimize(history_d, violation_d, feasible_d, log_interval);

	opt = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter, tol);
	opt->initialize_parallel();
	opt->optimize_parallel(history_d_p, violation_d_p, feasible_d_p, log_interval);

	std::unique_ptr<Optimizer<dimension>> opt1 = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter, tol, 0.9, 0.9, 1.2, 1.2, 0.3, 0.3);
	opt1->initialize();
	opt1->optimize(history_s, violation_s, feasible_s, log_interval);

	opt1 = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter, tol, 0.9, 0.9, 1.2, 1.2, 0.3, 0.3);
	opt1->initialize_parallel();
	opt1->optimize_parallel(history_s_p, violation_s_p, feasible_s_p, log_interval);

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
	*/
	return 0;
}

int serial_parallel_test(const typename OptimizerFactory<dimension>::OptimizerName &algorithm_name)
{
	// Initialize problem and solver parameters
	constexpr int log_interval = 50;
	int iter = 6000;
	int particles = 300;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	std::string alg_name_str = OptimizerFactory<dimension>::get_string_name(algorithm_name);

	// Preliminary informations to std out
	std::cout << "Serial vs Parallel optimization test for " << alg_name_str << std::endl;
	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num bees: " << particles << std::endl;

	// Test the serial version
	std::cout << std::endl;
	std::cout << "--- Serial optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<dimension>> opt = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter);
	auto t1 = std::chrono::high_resolution_clock::now();
	opt->initialize();
	opt->optimize();
	auto t2 = std::chrono::high_resolution_clock::now();
	opt->print_results();

	std::cout << "Exact value: " << TestProblems::get_exact_value<dimension>(test_problem) << std::endl;
	std::cout << std::setprecision(20) << "Absolute error: " << std::abs(TestProblems::get_exact_value<dimension>(test_problem) - opt->get_global_best_value()) << std::endl;
	std::cout << std::setprecision(20) << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - opt->get_global_best_position()).norm() << std::endl;
	double time_serial = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << std::setprecision(20) << "Elapsed initialization + optimization time: " << time_serial << std::endl;

	// Test the parallel version
	std::cout << std::endl;
	std::cout << "--- Parallel optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<dimension>> opt_p = OptimizerFactory<dimension>::create(algorithm_name, problem, particles, iter);
	t1 = std::chrono::high_resolution_clock::now();
	opt->initialize_parallel();
	opt->optimize_parallel();
	t2 = std::chrono::high_resolution_clock::now();
	opt->print_results();
	
	std::cout << "Exact value: " << TestProblems::get_exact_value<dimension>(test_problem) << std::endl;
	std::cout << std::setprecision(20) << "Absolute error: " << std::abs(TestProblems::get_exact_value<dimension>(test_problem) - opt->get_global_best_value()) << std::endl;
	std::cout << std::setprecision(20) << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - opt->get_global_best_position()).norm() << std::endl;
	double time_parallel = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << std::setprecision(20) << "Elapsed inizialization + optimization time: " << time_parallel << std::endl;
	std::cout << "Speedup: " << time_serial / time_parallel << std::endl;

	return 0;
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int mpi_rank = 0;
	int mpi_size = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (mpi_rank == 0)
	{
		std::cout << "Running with " << mpi_size << " MPI processes" << std::endl;
	}

	// Check the number of arguments
	if (argc != 3)
	{
		if (mpi_rank == 0)
		{
			std::cout << "Usage: ./test test_name algorithm_name" << std::endl;
			std::cout << "Available algorithms: ABC, SASPSO" << std::endl;
			std::cout << "Available tests for SASPSO algorithm: optimize static_adaptive, serial_parallel, time_numparticles" << std::endl;
			std::cout << "Available tests for ABC algorithm: optimize serial_parallel, time_numparticles" << std::endl;
		}
		MPI_Finalize();
		return -1;
	}
	// Create if it not exist the output directory
	if (mpi_rank == 0)
	{
		fs::create_directory("../output");
	}

	// Get from command line the algorithm to be tested
	std::string algorithm = argv[2];
	if (algorithm == "ABC")
	{
		// Get from command line the required test
		std::string test = argv[1];
		if (test == "optimize")
			optimize(OptimizerFactory<dimension>::ArtificialBeeColony);
		else if (test == "serial_parallel")
			serial_parallel_test(OptimizerFactory<dimension>::ArtificialBeeColony);
		else if (test == "time_numparticles")
			time_numparticles_test(OptimizerFactory<dimension>::ArtificialBeeColony);
		else
		{
			if (mpi_rank == 0)
			{
				std::cout << "Usage: ./test test_name algorithm_name" << std::endl;
				std::cout << "Available tests for ABC algorithm: optimize, serial_parallel, time_numparticles" << std::endl;
			}
		}
	}
	else if (algorithm == "SASPSO")
	{
		// Get from command line the required test
		std::string test = argv[1];
		if (test == "static_adaptive")
			static_adaptive_test(OptimizerFactory<dimension>::SelfAdaptiveSPSO);
		else if (test == "serial_parallel")
			serial_parallel_test(OptimizerFactory<dimension>::SelfAdaptiveSPSO);
		else if (test == "time_numparticles")
			time_numparticles_test(OptimizerFactory<dimension>::SelfAdaptiveSPSO);
		else if (test == "optimize")
			optimize(OptimizerFactory<dimension>::SelfAdaptiveSPSO);
		else
		{
			if(mpi_rank == 0)
			{
				std::cout << "Usage: ./test test_name problem_name" << std::endl;
				std::cout << "Available tests for SASPSO algorithm: optimize, static_adaptive, serial_parallel, time_numparticles" << std::endl;
			}
		}
	}
	else
	{
		if (mpi_rank == 0)
		{
			std::cout << "Usage: ./test test_name algorithm_name" << std::endl;
			std::cout << "Available algorithms: ABC, SASPSO" << std::endl;
		}
	}

	MPI_Finalize();

	return 0;
}