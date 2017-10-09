# VRP
Adaptive Heuristic Method Based on SA and LNS for Solving Vehicle Routing Problem

This program offers a heuristic solution to the Vehicle Routng Problem problem. The algorithm is based on two steps.
The first step uses an improved adaptive Simulated Annealing meta-heuristic to minimize the number of routes. 
The basic SA is modified to adaptively adjust its search span based on the quality of previous solutions obtained, 
i.e. if significant improvements do not occur repetitively, search span and depth are increased (and vice versa).
This adaptive feature significant improves the speed of algorithm.
The second part, inport the solution found via SA and minimizes the total transportation cost via a Large Neighborhood Search heuristic. 
The algorith is designed to handle very large distribution networks with multiple stores and dates. 
The algorithm also offers a simple solution to the shortcoming of many VRP solutions in optimizing fleet size. 
By adding a mechanism to assign routes to trucks the algorithm is capable of optimizing the fleet size as part of its objective function. 
