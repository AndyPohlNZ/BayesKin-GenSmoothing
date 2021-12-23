# BayesKin-GenSmoothing
This repository contains code and data to replicate the results described in the manuscript:
> Pohl A., Schofield M. & Ferber R. (2022). A generalised smoothing approach for continuous, planar, inverse kinematics problems

## Requirements 
R version 3.6.3 or higher and stan version 2.23 is required to run the scripts.  These are freely avaiable and we refer users to the relevant documentation for instalation at:

R: https://www.r-project.org/
Stan: https://mc-stan.org/

Users may wish to install the R studio IDE which is freely avaliable from: https://www.rstudio.com/

The models and analysis code presented rely on the following user contributed R pacakges:
- cmdstanr: https://mc-stan.org/cmdstanr/ 
- desolve: https://cran.r-project.org/web/packages/deSolve/index.html 
- splines2: https://cran.r-project.org/web/packages/splines2/index.html 
- numDeriv: https://cran.r-project.org/web/packages/numDeriv/index.html
- parallel: https://cran.r-project.org/web/views/HighPerformanceComputing.html 
- 
In addition the following packages are used for visulisation:
- latex2exp: https://cran.r-project.org/web/packages/latex2exp/index.html


Packages can be installed from the Comprehensive R archive Network (CRAN) by running the following in an R terminal:

    # Required Packages
    install.packages('desolve')
    install.packages('splines2')
    install.packages('numDeriv')
    install.packages('parallel')
    # Optional packages for visulaisation
    install.packages('latex2exp')

The `cmdstanr` package requires a sepearte installation process and this is outlined at https://mc-stan.org/cmdstanr/  


## Walkthrough
To run simulations open either `./gait-sim` to run the running gait examples or `./pend-sim` to run the pendulum examples.
Within each directory there is a `run.R` script which when run will generalte and solve 100 unique simulations.  If only one simulation is required modify set `nSim = 1` within the `run.R` script.

Once simulations have been completed, code within `analysis.R` generates the figures and results referenced in the main manuscript.

Further hints can be found within comments specific to each file.

## Directory contents

### pend-sim
- `library.R` Contains functions used within the main analysis.
- `run.R` Fits models to 100 sets of observed marker positions simulated from the
       pendulum model.
- `analysis.R` Generates results and plots found within the manuscript.

- `./model-files/`
    `sim-pendulum.stan` contains stan code to simulate the pendulum system
          given initial conditions
    `obs-pendulum.stan` contains stan code to simulate observed markers in the
           presence of measurement nosie sigma.

- `./model-fits/`
     contains R data files of each model fit to each simulation.

### gait-sim
- `library.R` Contains functions used within the main analysis.
- `run.R` Fits models to 100 sets of simulated marker positions
- `analysis.R` Generates results and plots found within the manuscript.

-  `./data/`
     `processes_data.R` processes raw marker tradjectories and static trial to a
         format suitiable of existing IK functions.
     `processed_data.rdata` contains the processed data from the motion capture
         experiement this is generated via the process_data.R script
     `s2_demographics.csv`, `s2_dynamicTrial.csv`, `s2_staticTrial.csv`, `s2_events.csv`
         contain the demographics, static and dynamic trial data from the
         motion capture experiement

- `./model-files/`
     `mod1.stan` contains stan code to fit the time independent BIK model
     `mod2.stan` contains stan code to fit the time continuous generalised smoothing BIK model

- `./model-fits/`
     contains R data files of each model fit to each simulation postfixed by
         the corresponding iteration's seed.

- `./gap-filling/`
     contains code to demonstrate the ability of the generalised smoothing approach
     to estimate kinematics in the presence of missing data




