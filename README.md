# UBOPF Case Studies
This repository contains a set of simulation case studies for the **Unbalanced Optimal Power Flow (UBOPF)** problem, applied to a **realistic low-voltage distribution network** with high penetration of **distributed photovoltaic (PV) generation**.

## Test Network Description

All scenarios are based on a common test system designed to represent a typical European-style residential distribution feeder. The main features of the network are:

- **18 nodes**
- **400 V line-to-line voltage**
- **Radial topology**
- **17 distribution lines** with typical parameters for low-voltage grids (high R/X ratio)
- **Unbalanced single-phase loads**, distributed across the three phases
- **Distributed photovoltaic generation** connected at various nodes and phases
- Modeled using a **three-phase unbalanced optimal power flow (UBOPF)** formulation implemented in Julia with JuMP.jl and Ipopt

This test network allows for realistic simulation of voltage rise, phase imbalance and inverter-based control strategies under diverse operating conditions.

## Folder Structure

All simulation cases are located in the `Case_studies/` folder, which includes:

- `Baseline_scenario/`  
  > Reference case with **low demand** and **high PV generation**.  
  > This scenario serves as a benchmark for evaluating voltage profiles and overvoltage risks.  
  > It also includes the analysis of two control strategies applied to the same network conditions.

- `Scenario_1/` to `Scenario_5/`  
  > Alternative operating scenarios designed to assess the impact of **three critical factors** on system performance:
  > 1. The **balance between demand and PV generation**  
  > 2. The **power factor of the loads**  
  > 3. The **degree of unbalance between phases**

### Script Structure

Each scenario folder (`Baseline_scenario/` to `Scenario_5/`) contains **three scripts**, corresponding to the control strategies evaluated:

1. `without_control.jl`
2. `P_control.jl`  
3. `PQ_control.jl` 

These scripts allow direct execution and comparison of results under the same network conditions and different control strategies.

## Control Strategies Evaluated

Across all case studies, the following PV inverter control strategies are analyzed:

- **No control strategy**  
- **P control**: Dynamic curtailment of active power 
- **PQ control**: Combined control of active and reactive power 

These strategies are integrated into the UBOPF formulation and tested for their ability to mitigate overvoltages and improve PV hosting capacity.

## Objectives

The main goals of these simulations are:

- To assess the occurrence and causes of **overvoltage phenomena** due to high PV penetration  
- To analyze the influence of key network factors on voltage profiles and control effectiveness  
- To compare the performance of different **inverter-based voltage control strategies**  
- To provide insights for the planning and operation of **unbalanced low-voltage networks** with high renewable integration


## Requirements

Make sure you have the following Julia packages installed:

```julia
using Pkg
Pkg.add("JuMP")
Pkg.add("Ipopt")
```

## Funding

This research was supported by the Madrid Government (Comunidad de Madrid-Spain) under the Multiannual Agreement 2023-2026 with Universidad Politécnica de Madrid, 'Line A - Emerging PIs' (grant number: 24-DWGG5L-33-SMHGZ1).
!(image.png)



