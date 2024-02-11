# SASPSO - Self Adaptive Standard Particle Swarm Optimization
The SASPSO algorithm is a modified version of the Standard PSO 2011 (SPSO 2011) algorithm that has been proposed by Tang et al. in [this](http://dx.doi.org/10.1155/2016/8627083) paper.

It is worthful to note that the SPSO 2011 algorithm, that has been proposed [here](https://ieeexplore.ieee.org/document/6557848) by Zambrano-Bigiarini et al., follows a similar idea
to the PSO algorithm implemented by us for the last assignment [here](https://github.com/AMSC22-23/PSO-marzo-santarsiero-guerrini) but it uses different strategies for some crucial steps.
For example the random point that has to be choosen in order to update the position is drawn uniformly inside a defined hypersphere, not an hypercube.

Another improvement is the dynamic adaptivity of the optimization parameters that the PSO algorithm needs. As shown by Tang et al., if we provide initial and final parameters that follows the given relations,
we can prove the convergence of the SASPSO 2011 algorithm.

This version can also solve contrained optimization problems (COPs) exploiting an adaptive relaxation method integrated with the feasibility-based rule. This technique is better than the common penalty rule since does not require additional parameters. More details can be found in [this](https://doi.org/10.1016/j.compchemeng.2011.09.018) paper by Zhang et al.

Two main classes have been developed:
- `Particle` represents a single particle of the swarm and manage its initialization and update at each iteration. It also implement the feasibility rule.
- `SASPSO` wraps the swarm of particles and the SASPSO 2011 algorithm besides with its initialization, execution, and output methods. 

The `SASPSO` class provides methods for the serial and the parallel implementation for both initialization and optimization methods.
The parallel implementation exploits OpenMP multithreading in order to speedup the particles update at each iteration, assigning to each thread a subset of the swarm.

Since there are no dependences between particles update at the same iteration, synchronization is only needed at the end of each iteration in order to implement s parallel reduction pattern
to find the best particle in the swarm.

## Required parameters
This algorithm requires only the basic parameters described by the common interface. Furter performance improvements can be achieved by tuning several optional parameters. In order to choose these
parameters the user should follow the indications in [this](http://dx.doi.org/10.1155/2016/8627083) paper.

## Available tests TODO
- `error_iteration`: Optimizes all the test functions with the same parameters both in serial and in parallel. Logs the error in function of the iterations count on the same optimizaiton loop for each function.Stores in the `error_iteration.csv` file all the errors in function of the iteration.
  > **Expected result**: It should show that the error decreases as the iteration count increases. The descent rate is random and starvation may be expected.

- `time_numparticles`: Optimizes several time the same test function with same parameters varying only the number of particles. The optimization is done both in serial and in parallel and the time is taken to analyse the speedup. Stores in the `time_numparticles.csv` file all the execution times of the serial and parallel optimize function for each swarm size of particles.
  > **Expected result**: It should show that the excution time increases linearly in both serial and parallel case. In particular we expect that the parallel increase rate is less than 
  the serial one. Teoretically the coefficient of the parallel case should be $\frac{1}{num\ threads}$ of the serial one. The speedup converges to the theoretical result as the number 
  of particles increases. It never reach that result due to the thread handling overhead.

- `serial_parallel_opt`: Basic optimization of a given test function both in serial and in parallel. Prints in standard output the execution time and the achieved error. Does not saves any file.
  > **Expected result**: It should show that the parallel version is faster than the serial one of about a factor of $num\ threads$. The error should be of the same magnitude, but may 
  vary due to the rendomness of the method.

## Documentation
The complete documentation of the public interface of our project can be consulted [here]() TODO.

## Author
Guerrini Alberto

## [<< Back to parent README](../../README.md)
