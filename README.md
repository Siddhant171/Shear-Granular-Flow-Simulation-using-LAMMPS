# Simulation-Based Analysis of Shear Stress in Granular Materials using LAMMPS

## Project Overview

This project investigates the shear stress behavior of dense granular assemblies as a function of shear rate using discrete element simulations in LAMMPS. The workflow generates physically consistent packings, compresses them to a dense state, and imposes shear at controlled rates to analyze resulting stress, velocity, and concentration fields.

The simulation combines discrete element modeling with MATLAB-based coarse-graining to interpret per-particle dynamics through continuum-level fields including normal/shear stresses, local shear rates, concentration gradients, and pressure distributions.

**Author**: Siddhantkumar Thute  
**Advisor**: Harkirat Singh  

## Problem Statement and Objectives

The objective is to characterize the shear stress behavior of a dense granular assembly as a function of shear rate. The central problem is understanding how the granular stress tensor and flow fields evolve across different shear rate cases.

Key objectives:
- Generate initial granular packing by pouring spherical particles under gravity
- Compress the poured configuration to a target pressure using wall-feedback algorithm
- Perform shear simulations at multiple shear rates (0.1, 0.2, 0.5 s⁻¹) and obtain per-particle kinematic and stress data
- Post-process data using MATLAB to compute bulk stresses, velocity gradients, concentration fields, and pressure profiles
- Compare results across shear rate cases and interpret stress response trends

## Methodology

### 1. Pouring Simulation (`in.pour`)
Creates initial granular packing with the following parameters:
- **Units**: SI (meters, kg, seconds, Pascals)
- **Geometry**: 2D simulation, sphere atom style
- **Boundaries**: Periodic in z, fixed in y, shrink-wrap in x
- **Particles**: Diameter d = 0.002 m, density 2500 kg/m³
- **Box Dimensions**: 0.1 × 0.032 m
- **Contact Model**: Granular Hooke with history
  - Normal stiffness: kₙ = 10⁵ N/m
  - Tangential stiffness: kₜ = 2×10⁴ N/m
  - Normal damping: γₙ = 50
  - Tangential damping: γₜ = 15
- **Dynamics**: NVE integration, gravity 9.81 m/s² downward, 2D enforcement
- **Output**: Restart file `restart3.poured` after 50,000 timesteps

### 2. Compression Stage (`in.compress`)
Compresses the poured assembly to target pressure of 100 Pa:
- **Input**: `Original_data.pour`
- **Contact Parameters**: Lower stiffness (kₙ = 400 N/m, kₜ = 0) for better packing
- **Pressure Calculation**: P_inst = -∑c_peratom[2] / (LₓLᵧ)
- **Feedback Control**: Proportional controller with gain 10⁻⁴, max displacement 5×10⁻⁶ m/cycle
- **Loop**: 10,000 iterations
- **Outputs**: 
  - `raw_pressure.dat`: Pressure time series
  - `dump_compress.dump`: Compression trajectory
  - `compressed_grains.data`: Final compressed configuration

### 3. Shear Simulation (`dump_after_eqm.shear`)
Performs shear deformation on compressed grains:
- **Shear Rates**: 0.1, 0.2, 0.5 s⁻¹ (in separate folders sr_0.1, sr_0.2, sr_0.5)
- **Geometry**: Triclinic box for simple shear
- **Boundaries**: Top/bottom grains fixed as shear plates, middle grains mobile
- **Contact Model**: Hooke with history (kₙ=400, kₜ=200, γₙ~1670, γₜ=50)
- **Pressure Control**: Maintains 100 Pa normal stress
- **Outputs per case**:
  - `dump_after_eqm.shear`: Equilibrium configuration
  - `raw_pressure.dat`: Pressure during shear
  - `sheared_compressed_grains.data`: Final sheared state
  - Atom dump files: Per-particle data (position, velocity, stress)

### 4. Data Extraction and Post-processing (`coarse_gran_2D.m`)
MATLAB script processes LAMMPS output to compute continuum fields:
- **Input Format**: Atom dumps with columns: id, type, x, y, z, vx, vy, vz, diameter, mass, c_peratom[1], c_peratom[2], c_peratom[4]
- **Bulk Stresses**:
  - Sₓₓ = (1/A) ∑ c_peratom[1]
  - Sᵧᵧ = (1/A) ∑ c_peratom[2]  
  - Sₓᵧ = (1/A) ∑ c_peratom[4]
- **Coarse-graining**: Lucy weighting kernel, window W=4d, N_slice=600
- **Computed Fields**:
  - Streamwise velocity Vₓ(z)
  - Shear rate ∇Vₓ(z)
  - Velocity curvature ∇²Vₓ(z)
  - Pressure P(z) = (Sₓₓ + Sᵧᵧ)/2
  - Pressure gradient ∇P(z)
- **Output**: stress_history_output.txt, visualization plots

### 5. Visualization
Three figures generated per shear rate:
1. Kinematics and concentration: Vₓ(z), γ̇(z)
2. Forces and scales: ∇²Vₓ(z), P(z), ∇P(z)
3. Stress evolution: Sₓₓ(t), Sᵧᵧ(t), Sₓᵧ(t), friction coefficient μ = Sₓᵧ/Sᵧᵧ

## Results and Analysis

### Shear Rate Dependence
The simulations show clear rate-dependent behavior:

- **Velocity Profiles**: Linear shear in bulk region, with magnitude scaling with imposed rate
- **Stress Response**: Bulk shear stress increases with shear rate
- **Pressure Profiles**: Sensitive to shear rate changes
- **Friction Coefficient**: Tracks transition to quasi-steady shearing regime

### Comparative Trends
- Increasing shear rate from 0.1 to 0.5 s⁻¹ scales bulk shear stress
- Local velocity profiles modify with rate
- Pressure gradients show rate sensitivity

## Files Description

### LAMMPS Scripts
- `in.pour`: Pouring simulation input
- `in.compress`: Compression with pressure feedback
- `dump_after_eqm.shear`: Shear simulation input

### MATLAB Analysis
- `coarse_gran_2D.m`: Post-processing and coarse-graining script

### Data Files
- `compressed_grains.data`: Compressed particle configuration
- `raw_pressure.dat`: Pressure time series
- `sheared_compressed_grains.data`: Sheared configuration
- `output-shear-*.txt`: Per-particle dump files
- `stress_history_output.txt`: Processed stress time series

### Documentation
- `report.pdf`: Complete project report with methodology, results, and analysis

## How to Run

1. **Pouring**: `lammps -in in.pour`
2. **Compression**: `lammps -in in.compress` 
3. **Shear**: `lammps -in dump_after_eqm.shear` (modify shear rate as needed)
4. **Analysis**: Run `coarse_gran_2D.m` in MATLAB

## Dependencies
- LAMMPS (with GRANULAR package)
- MATLAB
- OVITO (for visualization)

## Key Findings
- Captures rate-dependent granular rheology
- Extracts continuum fields from discrete simulations
- Demonstrates structured workflow for granular shear analysis
- Shows transition from solid-like to fluid-like behavior under shear

## Conclusion
The project successfully demonstrates a numerical workflow for granular shear simulation, capturing both local and bulk mechanical responses. The results indicate the pipeline can model rate-dependent stress and flow fields in dense granular assemblies.