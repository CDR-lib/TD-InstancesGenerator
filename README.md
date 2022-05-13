# Benchmarking Instances Generator for Aircraft Conflict Resolution
This project contains data files (.cpp code and user manual) for the generator presented in the paper M. Pelegrín & M. Cerulli (2021) "Aircraft conflict resolution: A benchmark generator".

Aircraft Conflict Resolution is one of the major tasks of computer-aided Air Traffic Management and represents a challenging optimization problem. This code generates benchmarking instances, each of them consisting of a set of flights with initial positions and vectors of velocities.

You will be able to generate:

	* Predefined 2D and 3D scenarios, including flow patterns and single encounter patterns.
	
	* Random 2D and 3D scenarios.
	
	* Pseudo-random 2D and 3D scenarios, where the trajectories meet a predefined traffic congestion.

The "src" folder contains generator's source (.cpp) code. The "scripts" folder contains bash scripts which can be executed to perform the computational study reported in Section 4 of the paper. Finally, in the "results" folder, several figures associated to different generable configurations are collected.
