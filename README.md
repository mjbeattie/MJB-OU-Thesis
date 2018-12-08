# MJB-Python-Git

This repository contains reusable Python code.  The files included do the following:

GA_beattie.py:  
This code is my modification of C. Nicholson's genetic algorithm code.
It solves a test function using GA heuristics.

PSO_beattie.py:  
This code is my codification of C. Nicholson's particle swarm
optimization code.  I have significantly changed the code to implement an
object-oriented approach.  This code allows the user to select the means by which
the next direction is chosen.  In particular, it contains the ability to use not
just past knowledge of the particle and swarm, but also the best "forward look".
The "forward look" is the gradient of the objective function seen by the particle
at a given iteration.

TABU_beattie.py:
This code implements a tabu search algorthim with best improvement.  It is used
for evaluating problems with discrete choice variable vectors, such as the
knapsack problem.

ANNEALING_beattie.py:
This code uses local search with simulated annealing to find an optimal solution.
This implementation solves the knapsack problem.

RANDOM_beattie.py:
This code implements the hill climbing with random walk algorithm.