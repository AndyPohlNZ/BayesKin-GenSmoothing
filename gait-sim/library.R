## # -------------------- Gait simulations--------------------------
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

# Load required packages
library("cmdstanr")
library("parallel")
library("splines2")
library("numDeriv")
library("latex2exp")

## --------------------------- Utilities ---------------
# Helper function to identify directory sturcture
stub <- function() {}
thisPath <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    if (length(grep("^-f$", cmdArgs)) > 0) {
        # R console option
        normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
    } else if (length(grep("^--file=", cmdArgs)) > 0) {
        # Rscript/R console option
        scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
    } else if (Sys.getenv("RSTUDIO") == "1") {
        # RStudio
        dirname(rstudioapi::getSourceEditorContext()$path)
    } else if (is.null(attr(stub, "srcref")) == FALSE) {
        # 'source'd via R console
        dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
    } else {
        stop("Cannot find file path")
    }
}

mod_pi <- function(angles) {
    # Rescales an angle on the (-inf, inf) scale to the interval [-180, 180)
    temp <- sapply(angles, function(ang) {
        if (abs(ang) <= pi) {
            return(ang)
        } else if (ang > pi) {
            remainder <- ang - pi
            return(-pi + remainder)
        } else if (ang < -pi) {
            remainder <- ang + pi
            return(pi + remainder)
        }
    })
    return(temp)
}
## Basis functions
fourier_basis <- function(x,
                          nK,
                          period,
                          intercept = TRUE,
                          derv = 0) {
    # Generates a fourier basis with nK terms with a period for the values in x.

    if (nK %% 2 != 0) {
        stop("nK must be even.")
    }
    omega <- (2 * pi) / period

    if (intercept) {
        basis <- matrix(0, nrow = length(x), ncol = nK + 1)

        ifelse(derv == 0,
            basis[, 1] <- 1,
            basis[, 1] <- 0
        )
        for (k in 1:(nK / 2)) {
            basis[, 2 * k] <- (k * omega)^derv *
                sin((k * omega * x) + (0.5 * derv * pi))
            basis[, 2 * k + 1] <- (k * omega)^derv *
                cos((k * omega * x) + (0.5 * derv * pi))
        }
    } else {
        basis <- matrix(0, nrow = length(x), ncol = nK)
        for (k in 1:(nK / 2)) {
            basis[, 2 * k - 1] <- (k * omega)^derv *
                sin((k * omega * x) + (0.5 * derv * pi))
            basis[, 2 * k] <- (k * omega)^derv *
                cos((k * omega * x) + (0.5 * derv * pi))
        }
    }
    return(basis)
}


bSpline_basis <- function(x,
                          nK,
                          degree = 3,
                          intercept = FALSE,
                          derv = 0) {
    # Generates a B-Spline basis of specified degree polynomial with nK+1 knots evenly spaced
    # throughout x.  if derv >1 then the derv deritive of the basis is returned.
    knots <- seq(min(x), max(x), length.out = nK + 1)[2:nK]
    basis <- bSpline(x, knots = knots, degree = degree, intercept = intercept)

    if (derv > 0) {
        basis <- dbs(x, knots = knots, degree = degree, intercept = intercept, derivs = derv)
    }

    return(basis)
}

# ---------------------- Sampling constants for MCMC----------------------------
# Set sampling contants
n_chains <<- 4 # number of MCMC chains
n_iter <<- 1000 # number of iterations after warmup
n_warmup <<- 1000 # number of warmup iterations
n_cores <<- 4
max_treedepth <<- 10

# ---------------------- Load pre-compiled Stan models -------------------
# compile (if necessary) stan models
mod1 <<- cmdstan_model(paste0(thisPath(), "/model-files/mod1.stan"), pedantic = TRUE)
mod2 <<- cmdstan_model(paste0(thisPath(), "/model-files/mod2.stan"), pedantic = TRUE)

## ---------------------- LS solution Functions --------------------------------------

gen_posture <- function(links,
                        r,
                        theta,
                        degrees = TRUE) {
    # Generates a posed system of links given translation of the origin (r) and
    # segment angles (theta)
    if (degrees) {
        theta <- theta * pi / 180
    }

    nlinks <- length(links)
    # rotation matrix
    Gam <- vector("list", nlinks)
    for (i in 1:nlinks) {
        Gam[[i]] <- matrix(c(cos(theta[i]), sin(theta[i]), -sin(theta[i]), cos(theta[i])), 2, 2)
    }

    # compute joint position
    J <- vector("list", nlinks + 1)
    J[[1]] <- matrix(r, 2, 1)
    for (i in 2:(nlinks + 1)) {
        J[[i]] <- J[[i - 1]] + Gam[[i - 1]] %*% c(links[[i - 1]]$length, 0)
    }

    # compute alpha
    alpha <- vector(mode = "list", length = nlinks)
    for (i in 1:nlinks) {
        nmarkers <- ncol(links[[i]]$x)
        alpha[[i]] <- matrix(NA, 2, nmarkers)
        for (j in 1:nmarkers) {
            x_ij <- links[[i]]$x[, j]
            alpha[[i]][, j] <- J[[i]] + Gam[[i]] %*% x_ij
        }
    }
    return(list(r = r, J = J, alpha = alpha))
}

cost <- function(params, y, links) {
    # Computes the least squares cost of a given pose given observations y
    # known links and parameters (params)
    # params = vector of parameters (rx, ry, theta1, theta2, ... thetan )
    nlinks <- length(links)
    r.hat <- params[1:2]
    theta.hat <- params[3:(nlinks + 2)]

    obs.hat <- gen_posture(links, r = r.hat, theta = theta.hat)
    y.hat <- obs.hat$alpha

    cost <- 0
    for (i in 1:nlinks) {
        nmarkers <- ncol(y.hat[[i]])
        for (j in 1:nmarkers) {
            diff <- y.hat[[i]][, j] - y[[i]][, j]
            cost <- cost + (sqrt(sum(diff * diff)))
        }
    }
    return(cost)
}


costgrad <- function(x, y, links) {
    # Numerical approximation of the LS cost function gradient at point x
    gradtemp <- grad(
        fun = cost,
        x = x,
        links = links,
        y = y,
        method.args = list(r = 6, zero.tol = 1e-6, eps = 1e-6, d = 1e-6)
    )
    return(gradtemp)
}

gen_LS_inits <- function(y,
                         links,
                         inits = c(0.0, 0.8, -90 * pi / 180, -90 * pi / 180, 0 * pi / 180)) {
    # Generates initial values via the time-indpenedent LS solution.

    # key parameters
    nL <<- length(links)
    nT <<- dim(y)[3]
    nM <<- dim(y)[2]
    r0 <- inits[1:2]
    theta0 <- inits[3:(nL + 2)]

    if (sum(is.na(y[, , 1, ])) > 0) {
        stop("y at t=0 cannot be missing")
    }

    # Compute first frame from inits
    r_inits <<- matrix(NA, nrow = nT, ncol = 2)
    theta_inits <<- matrix(NA, nrow = nT, ncol = nL)
    sigma_inits <<- rep(NA, nT)
    print(sprintf("Computing time: %.0f", 1))
    temp <- optim(
        par = c(r0, theta0 * 180 / pi), fn = cost, gr = costgrad,
        y = list(t(y[1, , 1, ]), t(y[2, , 1, ]), t(y[3, , 1, ])),
        links = links,
        method = "BFGS",
        control = list(trace = FALSE, maxit = 10000, reltol = 1e-8)
    )

    nmarkers <- 0
    obs.hat <- gen_posture(links, temp$par[1:2], temp$par[3:(nL + 2)])
    residual_sum <- 0
    for (i in 1:nL) {
        for (j in 1:nM) {
            residual_sum <- residual_sum +
                (y[i, j, 1, 1] - obs.hat$alpha[[i]][1, j])^2 + # xresidual
                (y[i, j, 1, 2] - obs.hat$alpha[[i]][2, j])^2 # yresidual
            nmarkers <- nmarkers + 1
        }
    }
    p <- length(temp)
    sigma_inits[1] <- sqrt(residual_sum / (2 * nmarkers))
    r_inits[1, ] <- temp$par[1:2]
    theta_inits[1, ] <- temp$par[3:(nL + 2)] * pi / 180

    # compute remaining frames iteratively using previous solution as new initial values
    for (t in 2:nT) {
        if (sum(is.na(y[, , t, ])) > 0) {

        }
        print(sprintf("Computing time: %.0f", t))
        temp <- optim(
            par = c(r_inits[t - 1, ], theta_inits[t - 1, ] * 180 / pi), fn = cost, gr = costgrad,
            y = list(t(y[1, , t, ]), t(y[2, , t, ]), t(y[3, , t, ])),
            links = links,
            method = "BFGS",
            control = list(trace = FALSE, maxit = 10000, reltol = 1e-8)
        )
        r_inits[t, ] <- temp$par[1:2]
        theta_inits[t, ] <- temp$par[3:(nL + 2)] * pi / 180

        nmarkers <- 0
        obs.hat <- gen_posture(links, temp$par[1:2], temp$par[3:(nL + 2)])
        residual_sum <- 0
        for (i in 1:nL) {
            for (j in 1:nM) {
                residual_sum <- residual_sum +
                    (y[i, j, t, 1] - obs.hat$alpha[[i]][1, j])^2 + # xresidual
                    (y[i, j, t, 2] - obs.hat$alpha[[i]][2, j])^2 # yresidual
                nmarkers <- nmarkers + 1
            }
        }
        p <- length(temp)
        sigma_inits[t] <- sqrt(residual_sum / (2 * nmarkers))
    }
    return(list(r_inits = r_inits, theta_inits = theta_inits, sigma_inits = sigma_inits))
}

## ---------------------- Simulate data function -------------------------------
simulate_data <- function(true_vals, links) {

    # Simulates marker observations given the 'true' segment angles
    # HJC location and noise
    # true_vals = list(r_true, theta_true, sigma_true)
    # links = link definitions from static trial.



    # Rotation matricies
    Gam <- array(NA, dim = c(nL, nT, nD, nD))
    for (t in 1:nT) {
        for (l in 1:nL) {
            Gam[l, t, 1:nD, 1:nD] <- rbind(
                c(cos(theta_true[t, l]), -sin(theta_true[t, l])),
                c(sin(theta_true[t, l]), cos(theta_true[t, l]))
            )
        }
    }

    # Joint centers
    J_sim <- array(NA, dim = c(nL, nT, nD))
    for (t in 1:nT) {
        J_sim[1, t, 1:nD] <- r_true[t, 1:nD]
        for (l in 2:nL) {
            J_sim[l, t, 1:nD] <- J_sim[l - 1, t, 1:nD] +
                (Gam[l - 1, t, 1:nD, 1:nD] %*% matrix(c(links[[l - 1]]$length, 0)))
        }
    }

    #  Marker Model
    alpha_sim <- array(NA, dim = c(nL, nM, nT, nD))
    y_sim <- array(NA, dim = c(nL, nM, nT, nD))
    for (t in 1:nT) {
        for (l in 1:nL) {
            for (m in 1:nM) {
                alpha_sim[l, m, t, ] <- J_sim[l, t, ] + Gam[l, t, , ] %*% x[l, m, ]
                y_sim[l, m, t, 1] <- rnorm(1, alpha_sim[l, m, t, 1], sigma_true[t])
                y_sim[l, m, t, 2] <- rnorm(1, alpha_sim[l, m, t, 2], sigma_true[t])
            }
        }
    }

    return(list(
        y_sim = y_sim, alpha_sim = alpha_sim, J_sim = J_sim
    ))
}


## ---------------------- RNG Functions --------------------------------------
rtruncnorm <- function(n = 1,
                       mean = 0,
                       sd = 1,
                       range = c(-1e6, 1e6)) {

    # Generates n random values from a truncated normal distribution (truncted
    # to 'range' via the inverse cdf method)
    x <-
        qnorm(runif(
            n,
            pnorm(range[1], mean = mean, sd = sd),
            pnorm(range[2], mean = mean, sd = sd)
        ),
        mean = mean,
        sd = sd
        )
    return(x)
}


rtrunccauchy <- function(n = 1,
                         mu = 0,
                         sigma = 1,
                         range = c(-1e6, 1e6)) {

    # Generates n random values from a truncated cauchy distribution (truncted
    # to 'range') via the inverse cdf method
    x <-
        qcauchy(runif(
            n,
            pcauchy(range[1], location = mu, scale = sigma),
            pcauchy(range[2], location = mu, scale = sigma)
        ),
        location = mu,
        scale = sigma
        )
    return(x)
}

## ---------------------- Initial value functions -----------------------------
# Sets initial values for MCMC sampling via LS solution.  Note that the ls_inits
# object must be in the local environment.
mod1_inits <- function() {
    # initial values for time-independent model - mod1.
    return(list(
        r = ls_inits$r_inits,
        theta = ls_inits$theta_inits,
        lsigma = log(ls_inits$sigma_inits)
    ))
}


mod2_inits <- function() {
    # Initial values for generalised smoothing model - mod2.
    rx_fit <- lm(ls_inits$r_inits[, 1] ~ bbasis - 1)
    ry_fit <- lm(ls_inits$r_inits[, 2] ~ fbasis - 1)
    theta1_fit <- lm(ls_inits$theta_inits[, 1] ~ fbasis - 1)
    theta2_fit <- lm(ls_inits$theta_inits[, 2] ~ fbasis - 1)
    theta3_fit <- lm(ls_inits$theta_inits[, 3] ~ fbasis - 1)
    beta_lsigma_fit <- lm(log(ls_inits$sigma_inits) ~ fbasis - 1)

    return(list(
        beta_rx = rx_fit$coef,
        mu_ry = ry_fit$coef[1],
        beta_ry = ry_fit$coef[2:(nK + 1)],
        mu_theta = c(
            theta1_fit$coef[1],
            theta2_fit$coef[1],
            theta3_fit$coef[1]
        ),
        beta_theta = rbind(
            theta1_fit$coef[2:(nK + 1)],
            theta2_fit$coef[2:(nK + 1)],
            theta3_fit$coef[2:(nK + 1)]
        ),
        delta = rnorm(nT, 0, 0.1),
        beta_lsigma = beta_lsigma_fit$coef,
        epsilon_lsigma = sqrt(mean((predict(beta_lsigma_fit) -
            log(ls_inits$sigma_inits))^2))
    ))
}

mod2_inits_missing_data <- function() {
    # Initial values for generalised smoothing model - mod2.  Modified to
    # handle missing data if gaps object is in the local environment.

    if (exists("gaps")) {
        midx <- c(gaps)
        rx_fit <- lm(ls_inits$r_inits[, 1] ~ bbasis[-midx, ] - 1)
        ry_fit <- lm(ls_inits$r_inits[, 2] ~ fbasis[-midx, ] - 1)
        theta1_fit <- lm(ls_inits$theta_inits[, 1] ~ fbasis[-midx, ] - 1)
        theta2_fit <- lm(ls_inits$theta_inits[, 2] ~ fbasis[-midx, ] - 1)
        theta3_fit <- lm(ls_inits$theta_inits[, 3] ~ fbasis[-midx, ] - 1)
        beta_lsigma_fit <- lm(log(ls_inits$sigma_inits) ~ fbasis[-midx, ] - 1)
        names(rx_fit$coef) <- ""
        names(ry_fit$coef) <- ""
        names(theta1_fit$coef) <- ""
        names(theta2_fit$coef) <- ""
        names(theta3_fit$coef) <- ""
        names(beta_lsigma_fit$coef) <- ""
    } else {
        rx_fit <- lm(ls_inits$r_inits[, 1] ~ bbasis - 1)
        ry_fit <- lm(ls_inits$r_inits[, 2] ~ fbasis - 1)
        theta1_fit <- lm(ls_inits$theta_inits[, 1] ~ fbasis - 1)
        theta2_fit <- lm(ls_inits$theta_inits[, 2] ~ fbasis - 1)
        theta3_fit <- lm(ls_inits$theta_inits[, 3] ~ fbasis - 1)
        beta_lsigma_fit <- lm(log(ls_inits$sigma_inits) ~ fbasis - 1)
    }
    return(list(
        beta_rx = rx_fit$coef,
        mu_ry = ry_fit$coef[1],
        beta_ry = ry_fit$coef[2:(nK + 1)],
        mu_theta = c(
            theta1_fit$coef[1],
            theta2_fit$coef[1],
            theta3_fit$coef[1]
        ),
        beta_theta = rbind(
            theta1_fit$coef[2:(nK + 1)],
            theta2_fit$coef[2:(nK + 1)],
            theta3_fit$coef[2:(nK + 1)]
        ),
        delta = rnorm(nT, 0, 0.1),
        beta_lsigma = beta_lsigma_fit$coef,
        epsilon_lsigma = sqrt(mean((predict(beta_lsigma_fit) -
            log(ls_inits$sigma_inits))^2))
    ))
}


#------------------------ Analysis Functions -----------------------------------
# The following functions aid in the generation of figures and numerical results
# referenced in the manuscript.

knee_joint_angle <- function(theta_thigh, theta_shank) {
    # Compute knee joint angle from thigh and shank angles.
    idx2 <- which(abs(theta_shank) > abs(theta_thigh))
    kja <- abs(theta_shank) - abs(theta_thigh)
    kja[idx2] <- -(abs(theta_thigh[idx2]) - abs(theta_shank[idx2]))
    return(kja)
}

ankle_joint_angle <- function(theta_shank, theta_foot) {
    # Compute ankle joint angle from shank and foot angles.
    idx2 <- which((theta_shank < 0) & (theta_foot > 0))

    aja <- abs(theta_foot) - abs(theta_shank) + (pi / 2)
    aja[idx2] <- pi - abs(theta_shank[idx2]) - abs(theta_foot[idx2]) - (pi / 2)

    return(aja)
}
posterior_mean <- function(samps_array) {
    # Compute the posterior mean given an array of MCMC samples
    if (length(dim(samps_array)) == 3) {
        dim(samps_array) <- c(dim(samps_array)[1] * dim(samps_array)[2], dim(samps_array)[3])
    }
    posterior_mean <- apply(samps_array, 2, mean)
    return(posterior_mean)
}

posterior_CI <- function(samps_array, quantiles = c(0.025, 0.975)) {
    # Compute the posterior credible interval given an array of MCMC samples
    # Set to compute central 95% credible interval by default
    if (length(dim(samps_array)) == 3) {
        dim(samps_array) <- c(dim(samps_array)[1] * dim(samps_array)[2], dim(samps_array)[3])
    }
    lwr <- apply(samps_array, 2, function(x) {
        return(quantile(x, probs = quantiles[1]))
    })
    upr <- apply(samps_array, 2, function(x) {
        return(quantile(x, probs = quantiles[2]))
    })
    return(cbind(lwr, upr))
}

## ---------------------- Model Performance -------------------------
# Evaluate model performance given simulated 'true' data r_sim, theta_sim and sigma_sim.

rmse <- function(y, yhat, minMax = TRUE) {
    # Root Mean Squared Error
    rmse <- sqrt(mean((y - yhat)^2))
    minRMSE <- min(sqrt((y - yhat)^2))
    maxRMSE <- max(sqrt((y - yhat)^2))

    if (minMax) {
        return(list(rmse = rmse, minRMSE = minRMSE, maxRMSE = maxRMSE))
    } else {
        return(rmse)
    }
}

compute_performance <- function(r_sim, theta_sim, sigma_sim, fit) {
    # computes performance metrics for a given model fit
    r1_draws <- fit$draws(paste0("r[", 1:nT, ",1]"))
    r2_draws <- fit$draws(paste0("r[", 1:nT, ",2]"))
    theta1_draws <- fit$draws(paste0("theta[", 1:nT, ",1]"))
    theta2_draws <- fit$draws(paste0("theta[", 1:nT, ",2]"))
    theta3_draws <- fit$draws(paste0("theta[", 1:nT, ",3]"))
    sigma_draws <- fit$draws(paste0("sigma[", 1:nT, "]"))
    nT <- nrow(r_sim)
    # RMSE
    rxrmse <- rmse(
        r_sim[, 1],
        posterior_mean(r1_draws)
    )
    ryrmse <- rmse(
        r_sim[, 2],
        posterior_mean(r2_draws)
    )
    theta1rmse <- rmse(
        theta_sim[, 1],
        posterior_mean(theta1_draws)
    )
    theta2rmse <- rmse(
        theta_sim[, 2],
        posterior_mean(theta2_draws)
    )
    theta3rmse <- rmse(
        theta_sim[, 3],
        posterior_mean(theta3_draws)
    )
    sigmarmse <- rmse(
        sigma_sim,
        posterior_mean(sigma_draws)
    )

    rmse_tab <- rbind(
        c(round(rxrmse$rmse * 1000, 3), round(rxrmse$minRMSE * 1000, 3), round(rxrmse$maxRMSE * 1000, 3)),
        c(round(ryrmse$rmse * 1000, 3), round(ryrmse$minRMSE * 1000, 3), round(ryrmse$maxRMSE * 1000, 3)),
        c(round(theta1rmse$rmse, 3), round(theta1rmse$minRMSE, 3), round(theta1rmse$maxRMSE, 3)),
        c(round(theta2rmse$rmse, 3), round(theta2rmse$minRMSE, 3), round(theta2rmse$maxRMSE, 3)),
        c(round(theta3rmse$rmse, 3), round(theta3rmse$minRMSE, 3), round(theta3rmse$maxRMSE, 3)),
        c(round(sigmarmse$rmse * 1000, 3), round(sigmarmse$minRMSE * 1000, 3), round(sigmarmse$maxRMSE * 1000, 3))
    )
    rownames(rmse_tab) <- c("rx[mm]", "ry[mm]", "theta1[rad]", "theta2[rad]", "theta3[rad]", "sigma[mm]")
    colnames(rmse_tab) <- c("mean", "min", "max")

    # CI Width
    rxCIWidth <- posterior_CI(r1_draws)[, 2] - posterior_CI(r1_draws)[, 1]
    ryCIWidth <- posterior_CI(r2_draws)[, 2] - posterior_CI(r2_draws)[, 1]
    theta1CIWidth <- posterior_CI(theta1_draws)[, 2] - posterior_CI(theta1_draws)[, 1]
    theta2CIWidth <- posterior_CI(theta2_draws)[, 2] - posterior_CI(theta2_draws)[, 1]
    theta3CIWidth <- posterior_CI(theta3_draws)[, 2] - posterior_CI(theta3_draws)[, 1]
    sigmaCIWidth <- posterior_CI(sigma_draws)[, 2] - posterior_CI(sigma_draws)[, 1]

    CI_tab <- round(cbind(
        c(
            mean(rxCIWidth) * 1000, mean(ryCIWidth) * 1000,
            mean(theta1CIWidth), mean(theta2CIWidth), mean(theta3CIWidth),
            mean(sigmaCIWidth) * 1000
        ),
        c(
            min(rxCIWidth) * 1000, min(ryCIWidth) * 1000,
            min(theta1CIWidth), min(theta2CIWidth), min(theta3CIWidth),
            min(sigmaCIWidth) * 1000
        ),
        c(
            max(rxCIWidth) * 1000, max(ryCIWidth) * 1000,
            max(theta1CIWidth), max(theta2CIWidth), max(theta3CIWidth),
            max(sigmaCIWidth) * 1000
        )
    ), 3)

    # Coverage
    rxCoverage <- mean((r_sim[, 1] > posterior_CI(r1_draws)[, 1]) & (r_sim[, 1] < posterior_CI(r1_draws)[, 2]))
    ryCoverage <- mean((r_sim[, 2] > posterior_CI(r2_draws)[, 1]) & (r_sim[, 2] < posterior_CI(r2_draws)[, 2]))
    theta1Coverage <- mean((theta_sim[, 1] > posterior_CI(theta1_draws)[, 1]) & (theta_sim[, 1] < posterior_CI(theta1_draws)[, 2]))
    theta2Coverage <- mean((theta_sim[, 2] > posterior_CI(theta2_draws)[, 1]) & (theta_sim[, 2] < posterior_CI(theta2_draws)[, 2]))
    theta3Coverage <- mean((theta_sim[, 3] > posterior_CI(theta3_draws)[, 1]) & (theta_sim[, 3] < posterior_CI(theta3_draws)[, 2]))
    sigmaCoverage <- mean((sigma_sim > posterior_CI(sigma_draws)[, 1]) & (sigma_sim < posterior_CI(sigma_draws)[, 2]))
    CI_tab <- cbind(CI_tab, round(c(rxCoverage, ryCoverage, theta1Coverage, theta2Coverage, theta3Coverage, sigmaCoverage) * 100, 1))

    rownames(CI_tab) <- c("rx[mm]", "ry[mm]", "theta1[rad]", "theta2[rad]", "theta3[rad]", "sigma[mm]")
    colnames(CI_tab) <- c("mean", "min", "max", "coverage")

    # Sampling efficiency
    rxsum <- fit$summary(paste0("r[", 1:nT, ",1]"))
    rysum <- fit$summary(paste0("r[", 1:nT, ",2]"))
    theta1sum <- fit$summary(paste0("theta[", 1:nT, ",1]"))
    theta2sum <- fit$summary(paste0("theta[", 1:nT, ",2]"))
    theta3sum <- fit$summary(paste0("theta[", 1:nT, ",3]"))
    sigmasum <- fit$summary("sigma")

    ESS_tab <- round(rbind(
        c(mean(rxsum$ess_bulk), min(rxsum$ess_bulk), max(rxsum$ess_bulk)),
        c(mean(rysum$ess_bulk), min(rysum$ess_bulk), max(rysum$ess_bulk)),
        c(mean(theta1sum$ess_bulk), min(theta1sum$ess_bulk), max(theta1sum$ess_bulk)),
        c(mean(theta2sum$ess_bulk), min(theta2sum$ess_bulk), max(theta2sum$ess_bulk)),
        c(mean(theta3sum$ess_bulk), min(theta3sum$ess_bulk), max(theta3sum$ess_bulk)),
        c(mean(sigmasum$ess_bulk), min(sigmasum$ess_bulk), max(sigmasum$ess_bulk))
    ), 0)
    rownames(ESS_tab) <- c("rx[n]", "ry[n]", "theta1[n]", "theta2[n]", "theta3[n]", "sigma[n]")
    colnames(ESS_tab) <- c("mean", "min", "max")
    ESS_sec_tab <- ESS_tab / (fit$time()$total)
    rownames(ESS_sec_tab) <- c("rx[n/s]", "ry[n/s]", "theta1[n/s]", "theta2[n/s]", "theta3[n/s]", "sigma[n/s]")
    colnames(ESS_sec_tab) <- c("mean", "min", "max")


    Rhat_tab <- round(rbind(
        c(mean(rxsum$rhat), min(rxsum$rhat), max(rxsum$rhat)),
        c(mean(rysum$rhat), min(rysum$rhat), max(rysum$rhat)),
        c(mean(theta1sum$rhat), min(theta1sum$rhat), max(theta1sum$rhat)),
        c(mean(theta2sum$rhat), min(theta2sum$rhat), max(theta2sum$rhat)),
        c(mean(theta3sum$rhat), min(theta3sum$rhat), max(theta3sum$rhat)),
        c(mean(sigmasum$rhat), min(sigmasum$rhat), max(sigmasum$rhat))
    ), 2)
    rownames(Rhat_tab) <- c("rx", "ry", "theta1", "theta2", "theta3", "sigma")
    colnames(Rhat_tab) <- c("mean", "min", "max")


    return(list(
        rmse = rmse_tab,
        interval = CI_tab, ess = ESS_tab,
        ess_sec = ESS_sec_tab,
        rhat = Rhat_tab,
        total_time = fit$time()$total
    ))
}

# ---------------------- Plotting Ultilities -------------------------
# Color plaette for plotting
myColors <<- list(
    grey = rgb(97 / 255, 97 / 255, 97 / 255, 0.2),
    orange = rgb(255 / 255, 127 / 255, 0 / 255, 1),
    blue = rgb(55 / 255, 78 / 255, 163 / 255, 1),
    red = rgb(228 / 255, 26 / 255, 28 / 255, 1),
    green = rgb(77 / 255, 175 / 255, 74 / 255, 1),
    purple = rgb(152 / 255, 78 / 255, 163 / 255, 1),
    pink = rgb(247 / 255, 129 / 255, 191 / 255, 1),
    violet = rgb(247 / 255, 129 / 255, 191 / 255, 1)
)

# Formatted plot labels
parmLabs <<- list(
    r = c(TeX("$r_x \\,  (m)$"), TeX("$r_y  \\, (m)$")),
    biasR = c(TeX("Bias: $r_x \\,  (m)$"), TeX("Bias: $r_y  \\, (m)$")),
    coverageR = c(TeX("Coverage: $r_x \\,  (%)$"), TeX("Coverage: $r_y  \\, (%)$")),
    intervalWidthR = c(
        TeX("95% credible interval width: $r_x \\, (mm)$"),
        TeX("95% credible interval width: $r_y \\, (mm)$")
    ),
    RMSER = c(
        TeX("$RMSE: $r_x (mm)$"),
        TeX("$RMSE: $r_y (mm)$")
    ),
    theta = c(
        TeX("$\\theta_1 \\, (rad))$"),
        TeX("$\\theta_2 \\, (rad)$"), TeX("$\\theta_3\\, (rad)) ")
    ),
    biasTheta = c(
        TeX("$Bias:\\, \\theta_1 \\, (rad)$"),
        TeX("$Bias:\\, \\theta_2 \\, (rad)$"),
        TeX("$Bias:\\, \\theta_3\\, (rad)$ ")
    ),
    RMSETheta = c(
        TeX("RMSE: $\\theta_1 \\, (rad)"),
        TeX("RMSE: $\\theta_2 \\, (rad)"),
        TeX("RMSE: $\\theta_3 \\, (rad)")
    ),
    coverageTheta = c(
        TeX("Coverage: $\\theta_1$ (%)"),
        TeX("Coverage: $\\theta_2$ (%)"),
        TeX("Coverage: $\\theta_3$ (%)")
    ),
    intervalWidthTheta = c(
        TeX("95% credible interval width: $\\theta_1 \\, (rad)$"),
        TeX("95% credible interval width: $\\theta_2 \\, (rad)$"),
        TeX("95% credible interval width: $\\theta_3 \\, (rad)$")
    ),
    sigma = TeX("$\\sigma \\, (m)$"),
    JA = c(TeX("Knee joint angle $(rad)$"), TeX("Ankle joint angle $(rad)$")),
    JA_bias = c(TeX("Bias: Knee joint angle $(rad)$"), TeX("Bias: Ankle joint angle $(rad)$")),
    JA_coverage = c(TeX("Coverage: Knee joint angle (%)"), TeX("Coverage: Ankle joint angle (%)")),
    JA_intervalWidth = c(TeX("95% credible interval width: Knee joint angle $(rad)$"), TeX("95% credible interval width: Ankle joint angle $(rad)$")),
    JA_RMSE = c(TeX("RMSE: Knee joint angle $(rad)$"), TeX("RMSE: Ankle joint angle $(rad)$"))
)