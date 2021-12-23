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
#
# ./figures/
#     contains figures generated for the manuscript

# --------------------- Perliminaries -------------------------------------
#  Load true data
load(file = "../data/processed_data.rdata")
y_true <- y

# Source functions in library file
source("../library.R")

# Number of simualtions to run

seed <<- 123 # Generate seeds for the simulation
set.seed(seed)
sfreq <<- 200 # Sampling frequency of simulated marker data

# ---------- Generate the LS solution as 'true - values' --------------------
ls_soln <- gen_LS_inits(y_true, links)
r_true <- ls_soln$r_inits
theta_true <- ls_soln$theta_inits
sigma_true <- ls_soln$sigma_inits

# Smooth true solution via lowess: Cleveland, W. S. (1979) Robust locally weighted regression and smoothing scatterplots. J. Amer. Statist. Assoc. 74, 829--836.
r_true[, 1] <- lowess(r_true[, 1] ~ ts, f = 0.1)$y
r_true[, 2] <- lowess(r_true[, 2] ~ ts, f = 0.1)$y
theta_true[, 1] <- lowess(theta_true[, 1] ~ ts, f = 0.1)$y
theta_true[, 2] <- lowess(theta_true[, 2] ~ ts, f = 0.1)$y
theta_true[, 3] <- lowess(theta_true[, 3] ~ ts, f = 0.1)$y
sigma_true <- lowess(sigma_true ~ ts, f = 0.1)$y


# ----------- Create the simulated observations including gaps ---------------
# Generate simulated data
sim_data <- simulate_data(
    true_vals = list(
        r_true = r_true,
        theta_true = theta_true,
        sigma_true = sigma_true,
    ),
    links = links
)
sim_data$ts <- ts


# create gaps
ngaps <- 1 # number of gaps to create
gap_length <- 20 * dT # 20 frames = 100ms
gaps <<- matrix(0, nrow = ngaps, ncol = gap_length / (dT)) # matrix of gap indicies
for (i in 1:ngaps) {
    g_t_range <- c((i - 1) * floor(nT / ngaps) + 1, i * floor(nT / ngaps) - gap_length / (dT))
    g_start <- sample(g_t_range[1]:g_t_range[2], size = 1)
    gaps[i, ] <- g_start:(g_start + (gap_length / (dT)) - 1)
}
# Remove gaps from simulated data
sim_data$y_sim <- sim_data$y_sim[, , -c(gaps), ]
sim_data$alpha_sim <- sim_data$alpha_sim[, , -c(gaps), ]
sim_data$J_sim <- sim_data$J_sim[, -c(gaps), ]
sim_data$ts <- sim_data$ts[-c(gaps)]

# Save true values
true_vals <- list(
    theta_true = theta_true,
    r_true = r_true,
    sigma_true = sigma_true,
    ts = ts,
    nT = nT,
    dT = dT,
    sampfreq = sampfreq,
    gaps = gaps,
    ngaps = ngaps
)
saveRDS(true_vals, "./true_vals.rdata")

# ----------------- Set up basis functions ----------------------------------
nK <<- 12
bbasis <<- bSpline_basis(ts, nK, degree = 3, intercept = TRUE, derv = 0)
fbasis <<- fourier_basis(ts, nK, period = range(ts)[2], intercept = TRUE, derv = 0)
dfbasis <<- fourier_basis(ts, nK, period = range(ts)[2], intercept = TRUE, derv = 0)


#------------------ Fit Generalised Smoothing BIK model ----------------------

# Generate LS solution for initial values.
print("    Generating LS initals")
ls_inits <<- gen_LS_inits(sim_data$y_sim, links)

# Process data into stan format
print("    Generating stan data for modelling")
N <<- nL * nM * length(sim_data$ts)
dat <- list()
dat$y <- matrix(NA, nrow = N, ncol = 2)
dat$ts <- dat$lidx <- dat$midx <- dat$tidx <- rep(NA, N) # initilise fields
idx <- 1
for (l in 1:nL) {
    for (m in 1:nM) {
        for (t in 1:length(sim_data$ts)) {
            dat$y[idx, 1:nD] <- sim_data$y_sim[l, m, t, ]
            dat$ts[idx] <- sim_data$ts[t]
            dat$lidx[idx] <- l
            dat$midx[idx] <- m
            dat$tidx[idx] <- which(ts == sim_data$ts[t])
            idx <- idx + 1
        }
    }
}

stan_data <- list(
    N = N,
    nL = nL,
    nM = nM,
    nD = nD,
    nT = length(ts),
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

# Fit model
print("    Fitting model 2")
mod2_fit <- mod2$sample(
    data = append(stan_data, list(
        nK = nK,
        nK2 = nK / 2,
        bbasis = bbasis,
        fbasis = fbasis,
        dfbasis = dfbasis
    )),
    init = mod2_inits_missing_data,
    chains = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    seed = seed,
    max_treedepth = max_treedepth,
    save_warmup = TRUE,
    parallel_chains = n_cores
)
print("        saving mod2")
mod2_fit$save_object(paste0("./mod2_gapfill_", seed, ".rds"))
print("        complete.")

# -------------- Generate Figures in manuscript -----------------------------------------------
# load fitted model and true value if not in environment
# mod2_fit <- readRDS("./mod2_gapfill_123.rds")
# true_vals <- readRDS("./true_vals.rdata")
# list2env(true_vals, .GlobalEnv)
# nT <<- nrow(r_true)
# nD <<- ncol(r_true)
# nL <<- ncol(theta_true)

# Compute performance metrics
mod2_perf <- compute_performance(r_true, theta_true, sigma_true, mod2_fit)

# Get draws from posterior for relevant parameters
rx_draws <- mod2_fit$draws(c(paste0("r[", 1:length(ts), ",1]")))
theta1_draws <- mod2_fit$draws(c(paste0("theta[", 1:length(ts), ",1]")))
theta2_draws <- mod2_fit$draws(c(paste0("theta[", 1:length(ts), ",2]")))
theta3_draws <- mod2_fit$draws(c(paste0("theta[", 1:length(ts), ",3]")))

# Determine knee and ankle joint angles.
kja <- aja <- array(NA, dim = dim(theta1_draws))
for (ii in 1:dim(theta1_draws)[1]) {
    for (jj in 1:dim(theta1_draws)[2]) {
        kja[ii, jj, ] <- knee_joint_angle(
            theta1_draws[ii, jj, ],
            theta2_draws[ii, jj, ]
        )
        aja[ii, jj, ] <- ankle_joint_angle(
            theta2_draws[ii, jj, ],
            theta3_draws[ii, jj, ]
        )
    }
}

# Plot knee and ankle joint angles with 95% credible intervals
pdf(
    width = 7 * 1, height = 1 * 7 * 0.4, pointsize = 12 * 0.5,
    file = paste0(thisPath(), "../figures/Gait_gapFill.pdf"),
    paper = "letter"
)
par(mfrow = c(1, 2))
plot(NA,
    xlim = range(ts), ylim = range(posterior_mean(kja)),
    xlab = "", ylab = "", bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
)
title(ylab = "Knee Joint Angle (rad)", line = 2.2, cex.lab = 1.5)
title(xlab = "Time (s)", line = 2.2, cex.lab = 1.5)
for (i in 1:ngaps) {
    rect(
        xleft = ts[min(gaps[i, ])], xright = ts[max(gaps[i, ])],
        ytop = pi,
        ybottom = -pi, col = "grey77", border = "grey88"
    )
}

lines(ts, posterior_mean(kja),
    col = "black", type = "l",
    bty = "n"
)
lines(ts, posterior_CI(kja)[, 1], col = "black", lty = 2)
lines(ts, posterior_CI(kja)[, 2], col = "black", lty = 2)
lines(ts, knee_joint_angle(true_vals$theta_true[, 1], true_vals$theta_true[, 2]),
    col = myColors$orange
)

plot(NA,
    xlim = range(ts), ylim = range(posterior_mean(aja)),
    xlab = "", ylab = "", bty = "n", , cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
)
title(ylab = "Ankle Joint Angle (rad)", line = 2.2, cex.lab = 1.5)
title(xlab = "Time (s)", line = 2.2, cex.lab = 1.5)
for (i in 1:ngaps) {
    rect(
        xleft = ts[min(gaps[i, ])], xright = ts[max(gaps[i, ])],
        ytop = pi,
        ybottom = -pi, col = "grey77", border = "grey88"
    )
}

lines(ts, posterior_mean(aja),
    col = "black", type = "l",
    bty = "n"
)
lines(ts, posterior_CI(aja)[, 1], col = "black", lty = 2)
lines(ts, posterior_CI(aja)[, 2], col = "black", lty = 2)
lines(ts, ankle_joint_angle(true_vals$theta_true[, 2], true_vals$theta_true[, 3]),
    col = myColors$orange
)
dev.off()