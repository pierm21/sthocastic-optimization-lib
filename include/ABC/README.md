### [<< Back to parent README](../../README.md)

# ABC - Artificial Bee colony   

The Artificial Bee Colony (ABC) algorithm is a metaheuristic optimization technique inspired by the foraging behavior of honey bees. It was introduced by Karaboga and Basturk and has proven effective for solving complex optimization problems, including constrained optimization problems (COPs).

The version developed in this project does refer to the Karaboga's ABC implementation for solving constrained optimization problem, which has been proposed [here](https://www.sciencedirect.com/science/article/pii/S1568494610003066).

The ABC algorithm simulates the intelligent foraging behavior of a honey bee swarm. The algorithm consists of three types of bees: employed bees, onlookers, and scouts, which work together to find the optimal solution.

## ABC key Features

- `Exploration and Exploitation`: The ABC algorithm balances exploration (global search) and exploitation (local search) through the cooperative behavior of the bees.
- `Adaptivity`: It dynamically adapts the search process based on the feedback from the environment, improving its ability to find optimal solutions.
- `Flexibility`: ABC can be applied to various types of optimization problems without significant modifications.

## Algorithm Description

It is worth to explain that, although different kind of bees are presented, they don't have to be considered effectively as if they are objects with different behaviours. In fact, in the algorithm, the employer phase, the scout phase and so on, just represent and characterize the different phases of the optimization single step. 

The ABC algorithm involves the following steps:

- `Initialization`: Generate an initial population of food sources (solutions) randomly.
- `Employed Bee phase`: Each employed bee evaluates the fitness of its food source and explores neighboring solutions, updating its position with a probability related to the algorithm parameter MR:

$$
x_{ij} = \begin{cases} 
x_{ij} + \beta_{ij}(x_{ij} - x_{kj}), & \text{if } R_j < MR \\
x_{ij}, & \text{otherwise}.
\end{cases}$$


- `Onlooker Bee Phase`: Onlooker bees select new food sources to be followed, based on their probability $p_i$, which is proportional to the fitness of the food source.

$$
p_i = \begin{cases} 
0.5 + \left( \frac{\text{fitness}_i}{\sum_{j=1}^{n} \text{fitness}_j} \right) \times 0.5, & \text{if solution is feasible} \\[12pt]
\left( \frac{1 - \text{violation}_i}{\sum_{j=1}^{n} \text{violation}_j} \right) \times 0.5. & \text{otherwise}
\end{cases}$$

- `Scout Bee Phase`: If a food source cannot be improved further, it is abandoned, and a scout bee randomly searches for a new food source, meaning that the bee is reinitialized.
- `Termination`: The process repeats until a stopping criterion (maximum number of iterations or acceptable solution quality) is met.

## Constrained Optimization

The ABC algorithm handles constrained optimization problems by incorporating a penalty function approach. The feasibility-based rule ensures that:
- feasible solutions are preferred over infeasible ones
- among feasible solutions, the one with better objective function value is chosen
- among infeasible solutions, the one with the smallest total constraint violation is preferred.

## Code structure

Two main classes have been developed:

-  `Bee` Represents a single bee in the swarm and manages its initialization, evaluation, and update at each iteration. It implements the feasibility rule and the search behavior of employed, onlooker, and scout bees.
-  `ABC` wraps the colony of bees and the ABC algorithm, providing methods for initialization, execution, and output. The class supports both serial and parallel implementations

The `ABC` class provides methods for the **serial** and the **parallel** implementation for both initialization and optimization methods. The parallel implementation uses OpenMP for multithreading and MPI for multiprocessing.
 
In the parallel version of the ABC algorithm, the entire colony of bees is divided equally among the available processors. This parallel implementation is designed for shared memory architectures and has been developed according to what proposed [here](https://ieeexplore.ieee.org/document/5393726), in fact some modifications with respect to the classical algorithm are needed in order to parallelize it. 

In particular, each thread handles a subset of the total bee's colony and operate completely on it, neglecting what is done by the other threads.

A final reduction, selecting the best bee among all the local colonies, will determine the global result.

TODO: MPI parallelism


## Tests and Results
Optimization of a Problem - `optimize`

This test optimizes a test function defined in the `test_problem` definte with the given dimensions. It logs the best value, total constraint violation, and the number of feasible solutions over iterations, storing the data in `output/abc_optimize_1.csv` file.

In the results below the G10 test problem in a 8D space is optimized. The swarm is composed by 5000 particles and 14k iterations are performed.
<p align="center">
  <img src="https://github.com/AMSC22-23/stochastic-optimization-lib/assets/48312863/f502d6f7-2b0f-4466-9eb8-18ac88ddd05b" height="500">
</p>

This plot highlights the policy to select the best among two particles that SASPSO 2011 utilizes:
1. A feasible solution is preferred over an infeasible solution
2. Among two feasible solutions, the one with better objective function value is preferred
3. Among two infeasible solutions, the one with smaller total constraint violation is chosen

Starting from a high-violation good-fitness solution and a very relaxed violation threshold the algorithm decreases the threshold converging to feasible solutions. The first non-feasible solution will have a better fitness, then when the 0 violation threshold is reached the first feasible solutions will have a very bad fitness. From now the constraint violation will be kept to 0, and the algorithms finds feasible solutions that minimizes the fiteness converging to the optimum.

Execution Time and Speedup - `time_numparticles`

This test logs the execution time of the ABC algorithm for varying swarm sizes, both in serial and parallel modes. It computes the parallel speedup and stores the data in output/abc_time_numparticles_num_threads.csv.

Conclusion

The ABC algorithm is a robust optimization method suitable for a wide range of optimization problems, including constrained ones. Its parallel implementation ensures scalability and efficiency, making it a valuable tool for solving complex problems.


### [<< Back to parent README](../../README.md)
