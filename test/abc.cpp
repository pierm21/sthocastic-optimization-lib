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

/*#define dimension 30
#define test_problem TestProblems::GRIEWANK
#define problem_name "GRIEWANK"*/

#define dimension 10
#define test_problem TestProblems::G7
#define problem_name "G7"


/*int optimize()
{
	constexpr int log_interval = 50;
	std::vector<double> best_values;
	int iter = 5000;
	int particles = 80;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	// Preliminary informations to std out
	std::cout << "Parallel adaptive optimization" << std::endl;
	std::cout << "Logs in /output/saspso_optimize.csv" << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Initialize history
	std::vector<double> history, violation, feasible, threshold;		

	// Get the exact global minimum
	double exact_value = TestProblems::get_exact_value<dimension>(test_problem);

	for (int i = 0; i < 1; i++)
	{
	std::unique_ptr<Optimizer<dimension>> opt_p = std::make_unique<ABC<dimension>>(problem, particles, iter);//, 1e-10);
	auto abc_ptr = dynamic_cast<ABC<dimension> *>(opt_p.get());
	if (!abc_ptr)
	{
		std::cerr << "Error: dynamic_cast failed" << std::endl;
		return 1;
	}

	//std::ofstream file_out;
	//file_out.open("../output/abc_optimize.csv");
	/*if (!file_out)
	{
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	// Write comments and header
	file_out << "# Fitness, constraints violation and feasible bees over iterations" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Bees: " << particles << std::endl;
	file_out << "# Problem: " << problem_name << std::endl;

	file_out << "Iters,value,violation,threshold,feasible_bees" << std::endl;

	//saspso_ptr->initialize_parallel();
	//saspso_ptr->optimize_parallel(history, violation, feasible, threshold, log_interval);
    abc_ptr->initialize_parallel();
	//abc_ptr->initialize();

    //abc_ptr->print_initizalization();
    abc_ptr->optimize_parallel();
	//abc_ptr->optimize();

    //abc_ptr->optimize(/*history, violation, feasible, threshold, log_interval);

    // Get the minimum value found
    double found_value = abc_ptr->get_global_best_value();
	/*RealVector<dimension> found_position =  abc_ptr->get_global_best_position();
	for (int i = 0; i < found_position.size(); i++)
	{
		std::cout << "Position[" << i << "]: " << found_position[i] << std::endl;
	}
	std::cout << "Exact value: " << exact_value << std::endl;
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

	//Print the final error
	//std::cout << std::endl << "Absolute error: " << std::abs(history.back() - exact_value) << std::endl;
	//std::cout << "Relative error: " << std::abs(history.back() - exact_value) / exact_value << std::endl;
    std::cout << std::endl << "Absolute error: " << std::abs(found_value - exact_value) << std::endl;
	if (exact_value != 0)
	{
	std::cout << "Relative error: " << std::abs(history.back() - exact_value) / exact_value << std::endl;
	}
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

	return 0;
}*/

int optimize()
{
	constexpr int log_interval = 50;

	int iter = 6000;
	int particles = 20;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	// Preliminary informations to std out
	std::cout << "Serial optimization" << std::endl;
	std::cout << "Logs in /output/abc_optimize.csv" << std::endl;

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num particles: " << particles << std::endl;

	// Initialize history
	std::vector<double> history, violation, feasible;

	std::unique_ptr<Optimizer<dimension>> opt_p = std::make_unique<ABC<dimension>>(problem, particles, iter);
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
	file_out << "# Fitness, constraints violation and feasible particles over iterations" << std::endl;
	file_out << "# Dimension: " << dimension << std::endl;
	file_out << "# Particles: " << particles << std::endl;
	file_out << "# Problem: " << problem_name << std::endl;

	file_out << "Iters,value,violation,feasible_particles" << std::endl;

	abc_ptr->initialize();
	abc_ptr->optimize(history, violation, feasible, log_interval);

	// Get the exact global minimum
	double exact_value = TestProblems::get_exact_value<dimension>(test_problem);

	// Print the final error
	std::cout << std::endl << "Absolute error: " << std::abs(history.back() - exact_value) << std::endl;
	if (exact_value != 0)
	{
	std::cout << "Relative error: " << std::abs(history.back() - exact_value) / exact_value << std::endl;
	}

	// Store on file
	for (int i = 0; i < history.size(); i++)
	{
		file_out << i * log_interval << ",";
		file_out << history[i] << ",";
		file_out << violation[i] << ",";
		file_out << feasible[i] << std::endl;
	}
	file_out.close();
	return 0;
}

// test the time as function of the number of particles
int time_numparticles_test()
{
	int iter = 6000;
	int max_particles = 1000;
	constexpr int log_interval = 20;
	constexpr int init_particles = 20;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	// Preliminary informations to std out
	std::cout << "Time and Speedup as function of colony size." << std::endl;
	std::cout << "Logs in /output/abc_time_numparticles.csv" << std::endl;

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
	file_out.open("../output/abc_time_numparticles.csv");
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
	std::cout << "Logging every " << log_interval << " iterations" << std::endl;

	for (int i = init_particles; i <= max_particles; i += log_interval)
	{
		// Print progress to stdout
	//	if ((i) % (log_interval) == 0)
			std::cout << "Starting test with " << i << " particle(s)" << std::endl;
		// Optimize the problem serially
		std::unique_ptr<ABC<dimension>> opt = std::make_unique<ABC<dimension>>(problem, i, iter);
		opt->initialize();
		auto t1 = std::chrono::high_resolution_clock::now();
		opt->optimize();
		auto t2 = std::chrono::high_resolution_clock::now();
		// Optimize parallel
		opt = std::make_unique<ABC<dimension>>(problem, i, iter);
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


int serial_parallel_test()
{
	// Initialize problem and solver parameters
	constexpr int log_interval = 50;
	std::vector<double> best_values;
	int iter = 6000;
	int particles = 80;
	auto problem = TestProblems::create_problem<dimension>(test_problem);

	std::cout << "Problem: " << problem_name << std::endl;
	std::cout << "Max iterations: " << iter << std::endl;
	std::cout << "Num bees: " << particles << std::endl;

	// Test the serial version
	std::cout << std::endl;
	std::cout << "--- Serial optimizer testing ---" << std::endl;
	std::unique_ptr<Optimizer<dimension>> opt = std::make_unique<ABC<dimension>>(problem, particles, iter);//, 1e-10);
	auto abc_ptr = dynamic_cast<ABC<dimension> *>(opt.get());
	if (!abc_ptr)
	{
		std::cerr << "Error: dynamic_cast failed" << std::endl;
		return 1;
	}	
	auto t1 = std::chrono::high_resolution_clock::now();
	abc_ptr->initialize();
	abc_ptr->optimize();
	auto t2 = std::chrono::high_resolution_clock::now();
	abc_ptr->print_results();

	std::cout << std::setprecision(20) << "Absolute error: " << std::abs(TestProblems::get_exact_value<dimension>(test_problem) - abc_ptr->get_global_best_value()) << std::endl;
	std::cout << std::setprecision(20) << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - abc_ptr->get_global_best_position()).norm() << std::endl;
	double time_serial = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << std::setprecision(20) << "Elapsed initialization + optimization time: " << time_serial << std::endl;

	// Test the parallel version
	std::cout << std::endl;
	std::cout << "--- Parallel optimizer testing ---" << std::endl;
std::unique_ptr<Optimizer<dimension>> opt_p = std::make_unique<ABC<dimension>>(problem, particles, iter);//, 1e-10);
	auto abc_p_ptr = dynamic_cast<ABC<dimension> *>(opt_p.get());
	t1 = std::chrono::high_resolution_clock::now();
	abc_p_ptr->initialize_parallel();
	abc_p_ptr->optimize_parallel();
	t2 = std::chrono::high_resolution_clock::now();
	abc_p_ptr->print_results();
	std::cout << std::setprecision(20) << "Absolute error: " << std::abs(TestProblems::get_exact_value<dimension>(test_problem) - abc_p_ptr->get_global_best_value()) << std::endl;
	std::cout << std::setprecision(20) << "Absolute distance: " << (TestProblems::get_exact_position<dimension>(test_problem) - abc_p_ptr->get_global_best_position()).norm() << std::endl;
	double time_parallel = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << std::setprecision(20) << "Elapsed inizialization + optimization time: " << time_parallel << std::endl;
	std::cout << "Speedup: " << time_serial / time_parallel << std::endl;

	return 0;
}

int main(int argc, char **argv)
{
	// Check the number of arguments
	if (argc != 2)
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for ABC algorithm: serial_parallel, time_numparticles" << std::endl;
		return -1;
	}
	// Create if it not exist the output directory
	fs::create_directory("../output");

	// Get from command line the required test
	std::string test = argv[1];
	if (test == "serial_parallel")
		serial_parallel_test();
	else if (test == "time_numparticles")
		time_numparticles_test();
	else if (test == "optimize")
		optimize();
	else
	{
		std::cout << "Usage: ./test-saspso test_name" << std::endl;
		std::cout << "Available tests for ABC algorithm: serial_parallel, time_numparticles" << std::endl;
		return -1;
	}

	return 0;
}