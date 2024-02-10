# Stochastic Optimization Library
Small collection of stochastic algorithms for constrained and unconstrained function optimization. The interface is common for each algorithm in order to facilitate the usage for the final user.


## Algorithms
Each of the implemented algorithms has its own ```README.md``` file in order to better explain its implementation, usage, implemented tests and their results.
- [**SASPSO - Self Adaptive Standard Particle Swarm Optimization**](include/SASPSO/README.md)
  - TODO: Short description
  - Constrained optimization available.
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

## Compile and Run
TODO: add informations

## Documentation
The full documentation is available [here](). TODO: add link to docs

## Authors
- Pierpaolo Marzo
- Domenico Santarsiero
- Alberto Guerrini
