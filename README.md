# Benchmarking Instances Generator for Aircraft Conflict Resolution
This project contains data files (.cpp code and user manual) for the generator presented in the paper M. Pelegrín & M. Cerulli (2021) "Aircraft conflict resolution: A benchmark generator".

Aircraft Conflict Resolution is one of the major tasks of computer-aided Air Traffic Management and represents a challenging optimization problem. This code generates benchmarking instances, each of them consisting of a set of flights with initial positions and vectors of velocities.

You will be able to generate:

	* Predefined 2D and 3D scenarios, including flow patterns and single encounter patterns.
	
	* Random 2D and 3D scenarios.
	
	* Pseudo-random 2D and 3D scenarios, where the trajectories meet a predefined traffic congestion.

The [src](https://github.com/CDR-lib/TD-InstancesGenerator/tree/main/src) folder contains generator's source (.cpp) code. The [scripts](https://github.com/CDR-lib/TD-InstancesGenerator/tree/main/scripts) folder contains bash scripts which can be executed to perform the computational study reported in Section 4 of the paper. Finally, in the [figures](https://github.com/CDR-lib/TD-InstancesGenerator/tree/main/figures) folder the figures associated to different generable configurations included in the paper are collected to illustrate the potential outputs of the generator.

# Generating instances and replicating results
When using the generator to produce instances to test, different parameters should be tuned. A detailed description of each of this parameters and their default values can be found in the [User Manual](https://github.com/CDR-lib/TD-InstancesGenerator/blob/main/USER-MANUAL.txt). 

If the generated instances are used as testbed, besides citing the paper "Aircraft conflict resolution: A benchmark generator", please specify to which values the different parameters are set. In this way, if anyone wants to use the same set of instances for testing purpose, it will only suffice to use the same values for the involved parameters.

# Citing 

# How to contribute
Thank you for considering contributing to our project! To report bugs and ask questions, please refer to the official [issue tracker](https://github.com/CDR-lib/TD-InstancesGenerator/issues). You can also address a problem by (1) forking the project, (2) correcting the bug / adding a feature, and (3) generating a pull request. However, we recommend that you first contact the authors and discuss the desired feature request.

You are also very welcome if you want to upload in our repository the instances you generated with our generator and used in your work! In this way, other researchers can use the same set of instances as benchmark, without having to re-generate such instances using the parameters you used. If this is the case, please send us an email and we will take care of this. 
