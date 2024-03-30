# Stochastic Optimization Library
Small collection of stochastic algorithms for constrained and unconstrained function optimization. The interface is common for each algorithm in order to facilitate the usage for the final user.


## Algorithms
Each of the implemented algorithms has its own ```README.md``` file in order to better explain its implementation, usage, implemented tests and their results.
- [**SASPSO - Self Adaptive Standard Particle Swarm Optimization**](include/SASPSO/README.md)
  - The SASPSO algorithm is a modified version of the Standard Particle Swarm Optimization 2011 (SPSO 2011) algorithm with **adaptive search parameters**.
  - Constrained optimization available with a constraint handling technique based on the **adaptive relaxation method** integrated with the **feasibility-based rule**.
- **SECOND ALGORITHM**

## Common interface
All the provided algorithms shares the same standard public interface. The user can use the same sequence of method calls to solve a given problem, indipendently from the choosen algorithm.

However, each algorithm gives the possibility to specify some additional parameters that are different from one algorithm to another given their different nature.
Custom parameters are preferred when solving specific and/or hard optimization problems in order to gain in performances.

The ```Problem``` class describes an optimization problem to be optimized. The user needs to instantiate a ```Problem``` object and set the space dimension in which it belongs, its fitness function, constraits, and bounds.

The ```Optimizer``` abstract class describes the interface for a specific optimizer implementation. The user needs to instantiate an object of the preferred algorithm passing the Problem as parameter.

The following code implements a possible basic usage of an optimizer.
```
TODO: sample code
```

## Requirements
- CMake
- C++17
- OpenMP 3.1
- Python3
- Eigen 3.3

NB: Eigen 3.3 can be loaded using [mk modules](https://github.com/pcafrica/mk) or downloading it from the official website.

## Compile and Run
1. Create the folders needed for building the project and saving the output `.csv` and `.png` files
   ```
   mkdir build
   mkdir output
   ```
2. Move into the build folder
   ```
   cd ./build
   ```
3. Execute cmake (it may require the flag with the explicit path to the Eigen library if it isn't in the default mk modules folder)
   ```
   cmake ..
   ```
4. Build the project
   ```
   cmake --build .
   ```
5. Launch the test executable for the preferred algorithm. A suitable test name must be provided (a complete list of the available test is shown if no parameter is provided)
   ```
   ./test-algorithm [test_name]
   ```
   The execution of a test may produce an output file `algorithm_test_name.csv` that can be found under the `output` folder.
6. Plot the results in a graphical way using the `csv_ploter.py` script passing as argument only the filename of the .csv file stored in the the `output` folder.
   ```
   cd ../scripts
   python csv_plotter.py algorithm_test_name.csv
   ```

## Documentation
The full documentation is available [here](). TODO: add link to docs

## Authors
- Pierpaolo Marzo
- Domenico Santarsiero
- Alberto Guerrini
