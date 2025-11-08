# mpi-battery-heat  

## Overview  

This project implements a **parallel 3D heat solver with Arrhenius chemical source term**, written in **C++17 with MPI**, and accompanied by **Python**.  

The solver models **thermal runaway in battery-like materials**, solving the heat equation with a nonlinear source term across a distributed-memory domain.  
The project was carried out as part of **scientific computing and HPC practice**, with emphasis on MPI parallelism, and numerical stability.  

## Goals  

- Implement **explicit (RK2)** and **implicit-explicit (IMEX)** time integration schemes.  
- Explore **stiff source terms** (Arrhenius-type kinetics) in PDE solvers.

## Approach  

- The project is structured with separation of concerns:  
  - `include/mbh/` → headers implementing grid, halo exchange, Laplacian, source terms, IO, and integrators.  
  - `src/` → main program logic.  
  - `scripts/` → Python plotting utilities (`plot_slice.py`, `plot_maxT.py`).  
  - `CMakeLists.txt` → build configuration.  
- Parallelism is handled via **MPI Cartesian topology** with halo exchanges.  
- Input handled through **CLI arguments** parsed into `Params`.  
- Outputs written as both **binary snapshots** and **CSV slices** for visualization.  

### Implementation Details  

- Default domain: **96 × 96 × 96 grid**, 64 mm³ cube.  
- Physics parameters:  
  - Density: 2500 kg/m³  
  - Heat capacity: 1000 J/kgK  
  - Conductivity: 1.0 W/mK  
  - Arrhenius parameters: A = 1.0e7, Ea = 8.0e4 J/mol, H = 3.0e5  
- Time integration schemes:  
  - **RK2** (explicit, fast but limited by CFL stability).  
  - **IMEX** (stiff-aware, stable with larger timesteps).  
- Initial condition: uniform background + **localized hotspot** (temperature, radius, position configurable). 

## Key Results  

- Correct scaling of explicit vs IMEX schemes.  
- Demonstrated **localized heating → spreading → runaway** under Arrhenius kinetics.  
- MPI parallelism tested on up to 8 ranks with domain decomposition.  
- Postprocessing tools successfully generated slice plots and tracked maximum temperatures.  

## Conclusions  

- The solver demonstrates the **trade-offs between explicit and implicit schemes** for stiff thermal problems.  
- Even a minimal MPI-based solver can handle **3D PDEs with nonlinear sources** efficiently.  
- The modular design (grid, Laplacian, source, integrators) supports **future extension** to more advanced PDEs.  
- The Python utilities enable straightforward visualization and diagnostics of HPC simulations.  

## How to Run  

### Requirements  

- **MPI** (OpenMPI or MPICH)  
- **CMake ≥ 3.16**  
- **C++17 compiler** (GCC, Clang)  
- **Python 3** with `numpy`, `matplotlib`  

### Build  

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

### Run Simulation  

```bash
# 64³ grid, RK2 scheme, 4 MPI ranks
mpirun -np 4 ./battery_heat --nx 64 --ny 64 --nz 64 --dt 5e-5 --t_final 0.002 --scheme rk2

# IMEX scheme, larger hotspot
mpirun -np 8 ./battery_heat --nx 128 --ny 128 --nz 128 --scheme imex --hotspot_T 500 --hotspot_radius 0.01
```

#### CLI Options  

| Flag | Description | Default |
|------|-------------|---------|
| `--nx --ny --nz` | Grid resolution | 96 96 96 |
| `--Lx --Ly --Lz` | Physical domain size [m] | 0.064 0.064 0.064 |
| `--rho` | Density [kg/m³] | 2500 |
| `--cp` | Heat capacity [J/kgK] | 1000 |
| `--k` | Conductivity [W/mK] | 1.0 |
| `--dt` | Timestep | 5e-5 |
| `--t_final` | Final time | 0.01 |
| `--scheme rk2|imex` | Time integrator | rk2 |
| `--T0` | Initial temperature [K] | 300 |
| `--T_env` | Ambient temperature [K] | 300 |
| `--hotspot_T` | Hotspot initial temperature [K] | 380 |
| `--hotspot_radius` | Hotspot radius [m] | 0.004 |
| `--hotspot_x --hotspot_y --hotspot_z` | Hotspot center (fraction of domain) | 0.5 0.5 0.5 |
| `--enable_source` | Enable Arrhenius source | true |
| `--A` | Arrhenius prefactor | 1.0e7 |
| `--Ea` | Activation energy [J/mol] | 8.0e4 |
| `--H` | Heat release per unit | 3.0e5 |
| `--outdir` | Output directory | out/ |
| `--px --py --pz` | MPI process grid | auto |

### Plotting

In main,

1. **2D slice visualization**  
   ```bash
   python scripts/plot_slice.py build/out/rank2_xy_mid80.csv --out slice.png
   ```

2. **Maximum temperature over time**  
   ```bash
   python scripts/plot_maxT.py "build/out/rank0_step*.bin"
   ```

---

## Acknowledgements  

HPC project by **Giorgos Kritopoulos**, studying parallel programming, numerical methods for PDEs. 

Date: 8 November 2025  
