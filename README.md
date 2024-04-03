# Stochastic Field Exploration of cluttered underwater Environments: Multi-Robot Informative Path Planning under Communication Constraints
This is a repo containing the code for my master's thesis. A paper focused on the source localization aspects of this work is currently under review for ICRA 2024.

## Content
This repo provides a comprehensive approach to the task of autonomous underwater field exploration using a team of multiple AUVs. For this, three independent libraries are implemented:
- **gmrf**: Contains an efficient implementation of [Joint estimation of gas and wind maps for fast-response applications](https://doi.org/10.1016/j.apm.2020.06.026) for underwater scenarios,
 extended to account for continuous measurement locations. Furthermore, an efficient source localization algorithm that leverages field estimates and is formulated as an MDP is implemented.
- **i_rrtstar**: Contains an efficient implementation of RRT* for the task of informative path planning. The planner can maximize the informativeness of the paths if a map of the expected information gain at different points is provided.
The planner optimizes the paths in a dedicated thread.
- **multiagent_sqp**: Contains a MultiAgent MPC problem that tracks the reference path, in this work provided by the RRT* planner, ensures continuous connectivity and collision avoidance between all agents, and maximizes individual information gain.
The problem is formulated as an SQP and then solved using CasaDi and OSQP.

Each library is configured using .toml files that can be found in the resources directory of the respective library.
## Dependencies
[SFML](https://www.sfml-dev.org/), [CasADi](https://web.casadi.org/) C++ (built with OSQP), and [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page).
[LibInterpolate](https://github.com/CD3/libInterpolate).

sfml:
```
sudo apt install libsfml-dev
```
eigen3:
```
sudo apt install libeigen3-dev
```

libinterpolate:

see [LibInterpolate](https://github.com/CD3/libInterpolate)

casadi:
```
git clone https://github.com/casadi/casadi.git
cd casadi
mkdir build
cd build
cmake .. (configure with osqp with ccmake)
```
