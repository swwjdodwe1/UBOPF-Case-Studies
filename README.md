# UBOPF Case Studies
This repository provides a Julia-based implementation of a program designed to solve the Unbalanced Optimal Power Flow (UBOPF) problem in distribution networks. The tool is developed to support studies on the impact of distributed photovoltaic (PV) generation in low-voltage systems.

# Repository Structure
1. Codigo_programa_UBOPF.jl
Main file that contains a general formulation of the UBOPF problem, implemented in the Julia programming language using the [JuMP.jl] package and the [Ipopt] solver.
This script applies the model to a simple 3-node test network, symbolically, to illustrate the basic operation of the program.

2. Simulaciones_UBOPF/ Folder
Contains the scripts corresponding to simulations carried out on a more realistic distribution network, composed of 18 nodes and distributed photovoltaic generation. This folder includes studies of different operating scenarios, as well as the evaluation of the impact of various control strategies.

Included Sub-Scenarios:
Baseline scenario: low demand and high photovoltaic generation.

Scenarios 1 to 5: variations in demand and PV generation profiles.

Evaluated Control Strategies:
No control strategies

Strategy 1: dynamic limitation of active power.

Strategy 2: reactive power control by photovoltaic inverters.

Each simulation produces results that allow analyzing the occurrence of overvoltages, the renewable generation hosting capacity, and the system’s behavior under different conditions.









Preguntar a ChatGPT
