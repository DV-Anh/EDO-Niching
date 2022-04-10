# EDO-Niching
The implementation of Niching Memetic Algorithm for Evolutionary Diversity Optimization on Traveling Salesperson Problem, as described in the work:<br>
Do, A.V., Guo, M., Neumann, A. and Neumann, F. (2022). Niching-based Evolutionary Diversity Optimization for the Traveling Salesperson Problem. Proceedings of the Genetic and Evolutionary Computation Conference. DOI: 10.1145/3512290.3528724

## How to use
Within the 2022_GECCO folder contains MATLAB implementation of NMA for EDO, as well as simple mutation-based EA for EDO. The chosen representation is visit-order permutation.

* **run.m**<br>The entry point where the experiment is run, containing settings with evaluation budgets, threshold values, etc.
* **tsp_instances.mat**<br>Contains data from 10 TSPLIB instances, including name, distance matrix, 2d Euclidean coordinates of vertices, and known optimal solution.
* **div_tsp_p1.m**<br>The NMA as a function.
* **dived.m**<br>(mu+1)-EA equalizing edge distances, maximizing sum-sum diversity.
* **divpd.m**<br>(mu+1)-EA maximizing smallest pairwise distances, maximizing sum-min diversity.

Run **run.m** as is to replicate the results. The output should be written into a separate file.
