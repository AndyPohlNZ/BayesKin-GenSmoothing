# # -------------------- Gait simulations--------------------------
# This collection of R scripts tests the performance of a time independent and generalised smoothing
# BIK models by evaluating their performance on marker observations simulated from
# a real marker based motion capture experiement

# Written as part of:
#      Pohl, Schofield and Ferber. (2022).
#      A generalised smoothing approach for continuous, planar, inverse kinematics problems
#
#       Code by: Andy Pohl
#       Human Perforamnce Labortory
#       Faculty of Kinesiology
#       University of Calgary
#       December 2021

# library.R Contains functions used within the main analysis.
# run.R Fits models to 100 sets of simulated marker positions
# analysis.R Generates results and plots found within the manuscript.

# ./data/
#     processes_data.R processes raw marker tradjectories and static trial to a
#         format suitiable of existing IK functions.
#     processed_data.rdata contains the processed data from the motion capture
#         experiement this is generated via the process_data.R script
#     s2_demographics.csv, s2_dynamicTrial.csv, s2_staticTrial.csv, s2_events
#         contain the demographics, static and dynamic trial data from the
#         motion capture experiement
#
# ./model-files/
#     mod1.stan contains stan code to fit the time independent BIK model
#     mod2.stan contains stan code to fit the time continuous generalised smoothing BIK model
#
# ./model-fits/
#     contains R data files of each model fit to each simulation postfixed by
#         the corresponding iteration's seed.
#
# ./gap-filling/
#     contains code to demonstrate the ability of the generalised smoothing approach
#     to estimate kinematics in the presence of missing data

# --------------------- Perliminaries -------------------------------------
#  Load true data
load(file = "./data/processed_data.rdata")
y_true <- y

# Source statements
source("./library.R")

nSim <<- 100 # Number of simualtions to run
seeds <<- 1:nSim # Generate seeds for each simulation
srate <<- 200 # Sampling rate of simulated marker data


# --------------------- Generate the LS solution to use as 'True Data' --------------------------
ls_soln <- gen_LS_inits(y_true, links)
r_true <- ls_soln$r_inits
theta_true <- ls_soln$theta_inits
sigma_true <- ls_soln$sigma_inits
# Downsample if srate is less than current samp frequency
if (srate > sampfreq) {
    stop("srate must be less than original sample frequency")
}
if (srate < sampfreq) {
    # Downsample ls_soln to given sample rate if necessary
    if (sampfreq %% srate != 0) {
        stop("srate must be a multiple of sampfreq")
    }
    downsampleFac <- sampfreq / srate
    downsample_idx <- seq(1, nT, by = downsampleFac)

    ts <- ts[downsample_idx]
    y_true <- y_true[, , downsample_idx, ]
    r_true <- r_true[downsample_idx, ]
    theta_true <- theta_true[downsample_idx, ]
    sigma_true <- sigma_true[downsample_idx]
    nT <<- length(ts)
    sampfreq <- srate
    dT <<- 1 / sampfreq
}

# Smooth true values via Lowess smoother
r_true[, 1] <- lowess(r_true[, 1] ~ ts, f = 0.1)$y
r_true[, 2] <- lowess(r_true[, 2] ~ ts, f = 0.1)$y
theta_true[, 1] <- lowess(theta_true[, 1] ~ ts, f = 0.1)$y
theta_true[, 2] <- lowess(theta_true[, 2] ~ ts, f = 0.1)$y
theta_true[, 3] <- lowess(theta_true[, 3] ~ ts, f = 0.1)$y
sigma_true <- lowess(sigma_true ~ ts, f = 0.1)$y

# Specify relevant parameters as a list
true_vals <- list(
    theta_true = theta_true,
    r_true = r_true,
    sigma_true = sigma_true,
    ts = ts,
    nT = nT,
    dT = dT,
    sampfreq = sampfreq
)
saveRDS(true_vals, "./true_vals.rdata") # Save true values to RDS file

# --------------------- Load true values if not in environent -----------------------------
# true_vals <- readRDS("./true_vals.rdata")
# list2env(true_vals, .GlobalEnv)
# nT <<- nrow(r_true)
# nD <<- ncol(r_true)
# nL <<- ncol(theta_true)

# ---------------------- Generate basis matricies for model-------------------------
# Basis matricies
nK <<- 12
bbasis <<- bSpline_basis(ts, nK, degree = 3, intercept = TRUE, derv = 0)
fbasis <<- fourier_basis(ts, nK, period = range(ts)[2], intercept = TRUE, derv = 0)
dfbasis <<- fourier_basis(ts, nK, period = range(ts)[2], intercept = TRUE, derv = 0)

# --------------------- Fit models----------------------------------------------
for (i in 1:nSim) {
    # Set seed for simulation
    seed <<- seeds[i]
    set.seed(seed)
    print(sprintf("Running iteration %.0f of %.0f", i, nSim))

    print("    Generating simulation")
    # Generate simulation
    sim_data <- simulate_data(
        true_vals = list(
            r_true = r_true,
            theta_true = theta_true,
            sigma_true = sigma_true,
        ),
        links = links
    )

    # Generate LS initials to initilize MCMC sampling
    print("    Generating LS initals")
    ls_inits <<- gen_LS_inits(sim_data$y_sim, links)

    # Convert data to stan format
    print("    Generating stan data for modelling")
    N <<- nL * nM * nT
    dat <- list()
    dat$y <- matrix(NA, nrow = N, ncol = 2)
    dat$ts <- dat$lidx <- dat$midx <- dat$tidx <- rep(NA, N)
    idx <- 1
    for (l in 1:nL) {
        for (m in 1:nM) {
            for (t in 1:nT) {
                dat$y[idx, 1:nD] <- sim_data$y_sim[l, m, t, ]
                dat$ts[idx] <- ts[t]
                dat$lidx[idx] <- l
                dat$midx[idx] <- m
                dat$tidx[idx] <- t
                idx <- idx + 1
            }
        }
    }
    stan_data <- list(
        N = N,
        nL = nL,
        nM = nM,
        nD = nD,
        nT = nT,
        ts = ts,
        lidx = dat$lidx,
        midx = dat$midx,
        tidx = dat$tidx,
        y = dat$y,
        x = x,
        ls = rbind(
            c(links$thigh$length, 0),
            c(links$shank$length, 0),
            c(links$foot$length, 0)
        )
    )

    # Fit Time independent BIK model
    print("    Fitting model 1")
    mod1_fit <- mod1$sample(
        data = stan_data,
        init = mod1_inits,
        chains = n_chains,
        iter_warmup = n_warmup,
        iter_sampling = n_iter,
        seed = seed,
        max_treedepth = max_treedepth,
        save_warmup = TRUE,
        parallel_chains = n_cores
    )
    print("        saving mod1")
    mod1_fit$save_object(paste0("./model-fits/mod1_", seed, ".rds"))
    print("        complete.")

    # Fit Generalised Smoothing BIK model
    print("    Fitting model 2")
    mod2_fit <- mod2$sample(
        data = append(stan_data, list(
            nK = nK,
            nK2 = nK / 2,
            bbasis = bbasis,
            fbasis = fbasis,
            dfbasis = dfbasis
        )),
        init = mod2_inits,
        chains = n_chains,
        iter_warmup = n_warmup,
        iter_sampling = n_iter,
        seed = seed,
        max_treedepth = max_treedepth,
        save_warmup = TRUE,
        parallel_chains = n_cores
    )
    print("        saving mod2")
    mod2_fit$save_object(paste0("./model-fits/mod2_", seed, ".rds"))
    print("        complete.")
}