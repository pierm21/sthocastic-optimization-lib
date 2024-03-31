### [<< Back to parent README](../../README.md)

# SASPSO - Self Adaptive Standard Particle Swarm Optimization
The SASPSO algorithm is a modified version of the Standard PSO 2011 (SPSO 2011) algorithm that has been proposed by Tang et al. in [this](http://dx.doi.org/10.1155/2016/8627083) paper.

It is worthful to note that the SPSO 2011 algorithm, that has been studied [here](https://ieeexplore.ieee.org/document/6557848) by Zambrano-Bigiarini et al., follows a similar idea
to the PSO algorithm implemented by us for the last assignment [here](https://github.com/AMSC22-23/PSO-marzo-santarsiero-guerrini) but it uses different strategies for some crucial steps.

For example in SPSO 2011 the random point that has to be choosen in order to update the position is drawn uniformly in a ball, not an hypercube. Moreover the parameter adaptivity does not require a tedious selection of the search parameters that would be required in any static PSO algorithm.

The adaptivity of the optimization parameters provided by SASPSO 2011, as shown by Tang et al., is done in order to have a proof for convergence given the initial and final parameters choice according to:

$$\begin{cases}
0 < \phi_{1s} + \phi_{1f} < 6ω_f + 6\\
−1 < \omega_f <\omega_s <1\\
\phi_{1s} = \phi_{2f} > \phi_{1f} = \phi_{2s} > 0
\end{cases}$$

Where $\omega$ is the inertia weight, $\phi_1$ is the cognitive acceleration parameter, and $\phi_2$ the social acceleration parameter. For each of them we denote the initial and final value as $x_i$ and $x_f$.

NB: The default parameters for the algorithm implemented in this library have been choosen according to these realations.

This version can also solve contrained optimization problems (COPs) exploiting an adaptive relaxation method integrated with the feasibility-based rule.
This technique is better than the common penalty rule since does not require additional parameters. More details can be found in [this](https://doi.org/10.1016/j.compchemeng.2011.09.018) paper by Zhang et al.

Two main classes have been developed:
- `Particle` represents a single particle of the swarm and manage its initialization and update at each iteration. It also implement the feasibility rule.
- `SASPSO` wraps the swarm of particles and the SASPSO 2011 algorithm besides with its initialization, execution, and output methods. 

The `SASPSO` class provides methods for the **serial** and the **parallel** implementation for both initialization and optimization methods.
The parallel implementation exploits OpenMP multithreading in order to speedup the particles update at each iteration, assigning to each thread a subset of the swarm.

At each iteration of the parallel version two synchronization points are needed:
- The first to ensure that the fesible particles count has been performed entirely before accessing the reduced varible.
- The second after the computation of the feasibility threshold in order to give the same updated value to all particles.

The reduction for the global best particle has been done using a critical section instead of a custom tree reduction since it has an easier implementation, the number of threads is usually limited, and the probability that all the threads arrives at the same instant is quite low due to the previous work sharing construct that has not a barrier at exit.

## Required parameters
This algorithm requires only the basic parameters described by the common interface. Furter performance improvements can be achieved by tuning several optional parameters. In order to choose these
parameters the user should follow the indications in [this](http://dx.doi.org/10.1155/2016/8627083) paper.

## Performed tests
The implementation of the SASPSO 2011 is provided with some tests in order to analyse its correctness and performance. The test source file for this algorithm is `test/saspso.cpp`.
In this section the provided tests are described and the results are shown.

### Optimization of a problem - `optimize`
This test optimizes the test function setted in the `test_problem` define with the given dimension. It logs the global best value, the related total constraint violation, the violation threshold, and the number of feasible particles in function of the iterations count, storing the data in the `output/saspso_optimize.csv` file.

In the results below the G10 test problem in a 8D space is optimized. The swarm is composed by 5000 particles and 14k iterations are performed.
<p align="center">
  <img src="https://github.com/AMSC22-23/stochastic-optimization-lib/assets/48312863/f502d6f7-2b0f-4466-9eb8-18ac88ddd05b" height="500">
</p>

This plot highlights the policy to select the best among two particles that SASPSO 2011 utilizes:
1. A feasible solution is preferred over an infeasible solution
2. Among two feasible solutions, the one with better objective function value is preferred
3. Among two infeasible solutions, the one with smaller total constraint violation is chosen

Starting from a high-violation good-fitness solution and a very relaxed violation threshold the algorithm decreases the threshold converging to feasible solutions. The first non-feasible solution will have a better fitness, then when the 0 violation threshold is reached the first feasible solutions will have a very bad fitness. From now the constraint violation will be kept to 0, and the algorithms finds feasible solutions that minimizes the fiteness converging to the optimum.

### Static vs adaptive optimization - `static_adaptive`
This test optimizes a given problem problem using the SASPSO 2011 algorithm using adaptive parameters and static ones to show the benefits of this dynamic algorithm. Moreover the optimization is performed both serially and in paralled to ensure no performance losses of the second one. 

In the results below the Gomez-Levy test problem in a 2D space is optimized. The swarm is composed by 5000 particles and 14k iterations are performed.

<p align="center">
  <img src="https://github.com/AMSC22-23/stochastic-optimization-lib/assets/48312863/0a5ea6e3-35a2-497a-b646-4701c670e42a" height="500">
</p>

NB: Small differences can derive from the intrinsic randomness of this algorithm. Note also that the constraint violation converges in few iterations to zero since the serach space have bigger feasibility areas than the one used before.

This plot shows that the parameter adaptivity feature of this implementation provides a much faster convergence with respect to static parameters (i.e. 200 vs 600 iterations). Both the serial and parallel version of the adaptive algorithm have comparable convergence rate. The adaptivity does not affect the convergence to the feasible search area as shown by the second plot.

### Execution time and Speedup - `time_numparticles`
This test optimizes several time a given test function varying only the number of particles. The optimization is done both serially and in parallel logging the execution time in order to compute the parallel speedup. The test stores in the `output/saspso_time_numparticles.csv` file all the logged execution time as function of the swarm size.

<p align="center">
  <img src="https://github.com/AMSC22-23/stochastic-optimization-lib/assets/48312863/cafccff2-9dda-4cea-8fbc-23a11b0fd172" height="500">
</p>

This result shows a parallel speedup around $9.25\times$ running on an Intel Core i7-13700H machine (14 hardware cores, 20 logical threads). The limited speedup is due to the non trivial synchronization between OpenMP threads needed at each iteration. The scalability of this implementation is linear with respect to the number of available cores.

## Author
Guerrini Alberto

### [<< Back to parent README](../../README.md)
