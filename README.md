# Tracking Interfaces in Random Logistic Free-Boundary Diffusion Problems: A Random Level Set Method

Reference (open access): 
**DOI:** [10.1007/s40314-025-03107-z](https://doi.org/10.1007/s40314-025-03107-z)  
Egorova, V., Casabán, M.C., Company, R. et al. Tracking interfaces in a random logistic free-boundary diffusion problems: a random level set method. Comp. Appl. Math. 44, 148 (2025). https://doi.org/10.1007/s40314-025-03107-z

## Abstract of the Paper

Free-boundary diffusive logistic model finds applications in diverse fields associated with population dynamics. These processes often possess stochastic characteristics and involve parameters with uncertainties. This study focuses on enhancing a two-dimensional diffusive logistic partial differential model with free boundary by incorporating randomness in the mean square sense, considering the conditions for well-posedness in the random case, which is crucial for the further analysis. Both unknown stochastic processes the solution and its moving front, and the parameters involved in the random problem as random variables, are constrained by a finite degree of randomness. To tackle this challenge, we propose a random level set method. Given the complexity of the problem, we employ alternating direction explicit methods for the interior solvers, to effectively address computational challenges. Since computing the mean and the standard deviation of both unknown stochastic processes are required, we combine the sample approach of the difference schemes together with Monte Carlo technique avoiding the storage accumulation of symbolic expressions of all the previous levels of the iteration process. Parallel computing is employed to enhance performance. A careful numerical analysis is performed in the mean square context to ensure stability, positivity, and boundedness. The set of presented examples illustrates these qualitative properties, assess numerical convergence and enables us to gain a deeper understanding of the system’s behavior attending to the geometry of the initial habitat. This approach provides valuable tools for analyzing and predicting spreading-vanishing dichotomy.

## Key Features

*   **Random Level Set Method (RLSM):** Implements the RLSM to track evolving interfaces in a stochastic environment.
*   **Free-Boundary Problem:** Models population dynamics where the habitat boundary is not fixed but part of the solution.
*   **Stochastic Modeling:** Incorporates randomness into model parameters (e.g., diffusion coefficient `D` and front expansion rate `μ`) using a Monte Carlo approach.
*   **Multiple Geometries:** Supports various initial habitat shapes to study the impact of geometry on population spreading:
    *   Ellipse (and Circle as a special case)
    *   Square (and Rectangle as a special case)
    *   Equilateral Triangle
*   **Comparison of Solvers:** Implements and allows comparison of several explicit finite difference schemes:
    1.  **EFDM:** Classic Explicit Finite Difference Method.
    2.  **IEFDM:** Improved Explicit Finite Difference Method (Larkin, 1964).
    3.  **ADE:** Alternating Direction Explicit (Saul'ev) Method.
    4.  **RALADE:** Random Average Larkin Alternating Direction Explicit Method (most stable and accurate).
*   **Visualization:** Generates plots for population density, front evolution, and total population dynamics over time.

  ## Requirements

*   **MATLAB R2021a or newer.** The code was developed and tested on MATLAB R2024a.
*   (Optional) **Parallel Computing Toolbox** is recommended for running a large number of Monte Carlo simulations (`Nsim > 1`) to improve performance. The provided code uses a standard `for` loop, which can be converted to a `parfor` loop for parallel execution.

## How to Run the Simulation

1.  Clone this repository or download the source files.
2.  Open the main script `main_simulation.m` in MATLAB.
3.  **Configure the simulation parameters** at the beginning of the script:

    ```matlab
    % --- Main Configuration ---
    
    % 1: Ellipse, 2: Square, 3: Triangle, 4: Rectangle
    form = 1; 
    
    % 1: EFDM, 2: IEFDM, 3: ADE Sauliev, 4: RALADE (Recommended)
    method = 4;
    
    % Number of Monte Carlo realizations. Set to 1 for a deterministic run.
    Nsim = 100;
    
    % Set to 1 to save figures, 0 to disable.
    save_flag = 0; 
    
    % Other parameters (Total time T, grid size Nx, Ny, etc.)
    T = 3;
    Nx = 200; 
    ```

4.  **Run the script.** The code will execute the simulation and generate several figures showing:
    *   The initial population and habitat shape.
    *   The evolution of the habitat front (initial vs. final mean front).
    *   The population dynamics `P(t)` for each Monte Carlo run, along with the mean.
    *   The final mean population distribution.

### Important Notes

*   **Deterministic vs. Stochastic:** To run a single, deterministic simulation, set `Nsim = 1`.
*   **Random Data:** For stochastic simulations (`Nsim > 1`), the code is configured to load pre-generated random variables from `.mat` files (e.g., `Dw_K100_std001.mat`). If you wish to run with a different `Nsim` value, you will need to generate and save your own random variables using script `GenerateRandomVariables_D_mu.mlx`
*   **Saving Figures:** If `save_flag` is set to `1`, you must update the `fname` variable to a valid directory path on your local machine.

## Repository Structure

```
.
├── main_simulation.m                      # The main MATLAB script to run simulations
├── Dw_K100_std001.mat                     # Example data file with random variables for Nsim=100
├── GenerateRandomVariables_D_mu.mlx       # The script to generate random parameters
└── README.md                              # This README file
```

## Citation

If you use this work or code in your research, please cite our paper:

Egorova, V., Casabán, M.C., Company, R. et al. Tracking interfaces in a random logistic free-boundary diffusion problems: a random level set method. Comp. Appl. Math. 44, 148 (2025). https://doi.org/10.1007/s40314-025-03107-z
