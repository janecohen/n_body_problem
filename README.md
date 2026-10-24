# N-Body Gravitational Simulation

A Fortran program that simulates the gravitational motion of particles in two or three dimensions. It supports Leapfrog (kick–drift–kick), fourth-order Runge–Kutta (RK4), and eighth-order Runge–Kutta (RK8/DOP853 coefficients) integration, with optional backward integration.

## Requirements

- GNU Fortran (`gfortran`)
- Python with NumPy, SciPy, and Matplotlib **only** if regenerating RK8 coefficients or using the separate plotting scripts

## Compile and run

From the project directory, make a data directory and compile the Fortran program:

```bash
mkdir -p data
make
./main
```

When prompted, enter one number for each choice:

| Simulation | `1` = 2D N-body; `2` = 3D N-body; `3` = two-body circular orbit |
| Integrator | `1` = Leapfrog; `2` = RK4; `3` = RK8 |
| Direction | `1` = forward only; `2` = forward and backward |



## Change simulation parameters

Edit the corresponding `case` block in **`main.f90`**, then recompile. 

| Setting | Meaning |

| `N` | Number of particles (positive integer; use at least 2 for gravitational interactions) |
| `dim` | Spatial dimensions (`2` or `3`) |
| `iterations` | Number of integration steps (positive integer) |
| `dt` | Time step (positive real number) |
| `eps` | Gravitational softening parameter (nonnegative real number) |
| `omega` | Initial angular velocity (real number; use 1 as default) |
| `output_filename`, `quantity_filename` | Output text-file paths |
| `output_backwards_filename`, `quantity_backwards_filename` | Output paths for backward integration |

Use Fortran double-precision notation for real values, e.g. `dt = 0.001d0`. The output directories must already exist. 

## Output files

Results are written to the `data/` directory by default:

- `data_*.txt`: particle positions and velocities, one particle per row, for each saved time step.
- `quantity_*.txt`: potential energy, kinetic energy, angular momentum, and barycentric velocity for each saved time step.
- Files ending in `_backwards.txt`: corresponding to backward-integration results, when selected.

The program also prints the simulation settings and execution time in the terminal. Use the separate Python analysis scripts to visualize trajectories and assess conservation and computational cost.

## RK8 coefficients

The RK8 integrator reads `dop853_coefficients.txt` from the current working directory. Keep the supplied file next to the executable when running from the project directory. To regenerate it, run:

python3 rk8_coefficients.py
