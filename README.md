# Efficient_path_planning

**Overview:**
This repository implements a path-planning and task assignment algorithm for very large multi-robot systems that should fulfill a global Boolean specification.

There are 4 main programs as it follows:

1) main_map_reachability.m is for **reachability scenario** (number of regions to be reached is equal with the number of robots; for this simulations the environment, together with the initial and final positions of the robots are randomly generated);

Inputs: number of robots, number of regions;
Output: robot trajectories;

Note: 1) if cell capacity (s) should be greater than 1, in order to ensure collision avoidance the problem will be solved using (s) intermediary markings and cell capacity will be considered equal with 1.
2) The user can choose if plots with trajectories will be displayed.
3) The user can choose if the Integer Linear Programming (ILP) formulation will be also solved. This might take a while for large teams of mobile agents. 

2) main_map_benchmarks.m is for testing the algorithm on **4 benchmarks maps from [1]** (number of regions to be reached is equal with the number of robots; for this simulations the initial and final positions of the robots are randomly generated); the maps considered are room-32-32-4, random-32-32-20, den312d and ht_chantry;

Inputs: number of robots, number of regions;
Output: robot trajectories;

Note: 1) if cell capacity (s) should be greater than 1, in order to ensure collision avoidance the problem will be solved using (s) intermediary markings and cell capacity will be considered equal with 1.
2) The user can choose if plots with trajectories will be displayed.
3) The user can choose if the Integer Linear Programming (ILP) formulation will be also solved. This might take a while for large teams of mobile agents. 

3) main_map_boolean.m is for **global boolean-based goal** (number of regions to be reached is not necessarily equal with the number of robots; for this simulations the environment and the boolean formula, together with the initial and final positions of the robots are randomly generated);

Inputs: number of robots, number of regions;
Output: robot trajectories;

Note: 1) if cell capacity (s) should be greater than 1, in order to ensure collision avoidance the problem will be solved using (s) intermediary markings and cell capacity will be considered equal with 1.
2) The user can choose if plots with trajectories will be displayed.
3) The user can choose if the Integer Linear Programming (ILP) formulation will be also solved. This might take a while for large teams of mobile agents. 
4) It is not mandatory to have all the robots moving in order to fulfill the global specification.

3) main_map_smart_plant.m is for the **smart manufacturing plant scenario** (number of regions to be reached is not necessarily equal with the number of robots; for this simulations the environment and the boolean formula, together with the initial and final positions of the robots are randomly generated);

Inputs: number of robots, number of regions;
Output: robot trajectories;

Scenario: 50 robots evolving in a small smart plant manufacturing;
	Boolean specification defined over a number of 50 requests:
		- 10 jobs type 1 - inventory inspection and monitoring (should be 100% fulfilled) - the blue cells;
		- 15 jobs type 2 - hazard detection and mitigation (should be 100% fulfilled) - the cyan cells;
		- 10 jobs type 3 - product quality control (should be at least 50% fulfilled) - the dark green cells;
		- 15 jobs type 4 - maintenance and calibration (should be at least 75% fulfilled) - the yellow cells;

Simulation example:

![smart_plant_small_example](https://github.com/user-attachments/assets/d8edf394-9bb8-46ac-a21c-80ff8592303a)

Note: 1) if cell capacity (s) should be greater than 1, in order to ensure collision avoidance the problem will be solved using (s) intermediary markings and cell capacity will be considered equal with 1.
2) The user can choose if plots with trajectories will be displayed.
3) The user can choose if the Integer Linear Programming (ILP) formulation will be also solved. This might take a while for large teams of mobile agents. 
4) It is not mandatory to have all the robots moving in order to fulfill the global specification.

All problems are solved using intlinprog solver (for the relaxed LPs problems, intcon option is [ ], meaning that no variable is imposed to be integer);

**Note:** All simulations were performed on a computer equipped with an Intel Xeon Gold 6240 CPU at 2.60GHz and 1 TB RAM. Different perfomances may be obtained with different specifications.

[1] Stern, R., Sturtevant, N., Felner, A., Koenig, S., Ma, H., Walker, T., Li, J., Atzmon, D., Cohen, L., Kumar, T.K. and Barták, R., 2019. Multi-agent pathfinding: Definitions, variants, and benchmarks. In Proceedings of the International Symposium on Combinatorial Search (Vol. 10, No. 1, pp. 151-158).
