# # -------------------- Simple Pendulum Example --------------------------
# This example simulates observed marker positions from a simple pendulum
# model and implementes either a time independent, generalised smoothing
# or the true ODE model to generate solutions to the IK problem

# Written as part of:
#      Pohl, Schofield and Ferber. (2022).
#      A generalised smoothing approach for continuous, planar, inverse kinematics problems
#
#       Code by: Andy Pohl
#       Human Perforamnce Labortory
#       Faculty of Kinesiology
#       University of Calgary
#       December 2021

## library.R Contains functions used within the main analysis.
## run.R Fits models to 100 sets of observed marker positions simulated from the
#       pendulum model.
## analysis.R Generates results and plots found within the manuscript.

## ./model-files/
##     sim-pendulum.stan contains stan code to simulate the pendulum system
##          given initial conditions
##     obs-pendulum.stan contains stan code to simulate observed markers in the
#           presence of measurement nosie sigma.

## ./model-fits/
##     contains R data files of each model fit to each simulation.

# ---------------------- Import External Libraries -----------------------------
library("deSolve")
library("cmdstanr")
library("parallel")
library("splines2")
library("numDeriv")
library("latex2exp")

# ---------------------- Perliminaries---------------------------------
# Set sampling contants
n_chains <<- 4 # number of MCMC chains
n_iter <<- 1000 # number of iterations after warmup
n_warmup <<- 1000 # number of warmup iterations
n_cores <<- detectCores()

# ---------------------- Load pre-compiled Stan models -------------------
# compile (if necessary) stan models
sim_pend_mod <<- cmdstan_model("./model-files/sim-pendulum.stan", pedantic = TRUE)
gen_obs_mod <<- cmdstan_model("./model-files/obs-pendulum.stan", pedantic = TRUE)

mod1 <<- cmdstan_model("./model-files/mod1.stan", pedantic = TRUE)
mod2 <<- cmdstan_model("./model-files/mod2.stan", pedantic = TRUE)
mod3 <<- cmdstan_model("./model-files/mod3.stan", pedantic = TRUE)

# ---------------------- Utilities -----------------------
mod_pi <- function(angles) {
    # Rescales an angle on the (-inf, inf) scale to the interval [-pi, pi)
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

# ---------------------- Data Generation functions -----------------------
sim_pend <- function(ts, state0, parms) {
    # Simulate a simple pendulum over ts
    # given ICs state0 and parameters (g and length) in parms

    # Outputs:
    #   output: matrix[1:nT, 2] of theta and dtheta.

    nT <- length(ts)

    stan_data <- list(
        nT = nT - 1,
        state0 = state0,
        t0 = ts[1],
        ts = ts[-1],
        parms = c(parms$g, parms$Lm)
    )

    fit <-
        sim_pend_mod$sample(
            data = stan_data,
            chains = 1,
            iter_warmup = 0,
            iter_sampling = 1,
            seed = 1,
            save_warmup = FALSE,
            fixed_param = TRUE
        )
    fit

    theta_samps <- as.data.frame(fit$draws(paste0("state[", 1:(nT - 1), ",1]")))
    dtheta_samps <- as.data.frame(fit$draws(paste0("state[", 1:(nT - 1), ",2]")))
    tmp <- cbind(as.numeric(theta_samps), as.numeric(dtheta_samps[1, ]))
    output <- rbind(state0, tmp)
    return(output)
}


gen_link <- function(seg_length,
                     plate_placement,
                     marker_plate_dim = c(0.05, 0.1)) {
    # Creates a link of length seg.length with a marker plate at plate.center
    # with marker plate centerd at plate.placement.length*seg.length
    # Aligned on the positive x axis such that the origin of the link
    # is at c(0,0)

    plate_center <- plate_placement * seg_length
    plate_loc <- c(0, -plate_center)

    markers <-
        matrix(nrow = nD, ncol = nM)
    markers[1:2, 1] <- c(
        plate_loc[1] - marker_plate_dim[1] / 2,
        plate_loc[2] + marker_plate_dim[2] / 2
    )
    markers[1:2, 2] <- c(
        plate_loc[1] + marker_plate_dim[1] / 2,
        plate_loc[2] + marker_plate_dim[2] / 2
    )
    markers[1:2, 3] <- c(
        plate_loc[1] + marker_plate_dim[1] / 2,
        plate_loc[2] - marker_plate_dim[2] / 2
    )
    markers[1:2, 4] <- c(
        plate_loc[1] - marker_plate_dim[1] / 2,
        plate_loc[2] - marker_plate_dim[2] / 2
    )

    return(markers)
}

gen_obs <- function(ts, state0, x, r_true, sigma, parms) {
    # Solves the ODE for an simple pendulum at times ts
    # given initial states state0 and parms list
    #
    # Returns a list of marker observations y
    # and expected marker positions (alpha)

    nT <- length(ts)

    stan_data <- list(
        nT = nT,
        nM = nM,
        nD = nD,
        state0 = state0,
        t0 = ts[1],
        ts = ts,
        x = t(x),
        r = r_true,
        sigma = sigma,
        parms = c(parms$g, parms$Lm)
    )

    fit <- gen_obs_mod$sample(
        data = stan_data,
        chains = 1,
        iter_warmup = 0,
        iter_sampling = 1,
        seed = seed,
        save_warmup = FALSE,
        fixed_param = TRUE
    )

    y <- array(NA, dim = c(nM, nD, nT))
    alpha <- array(NA, dim = c(nM, nD, nT))
    for (m in 1:nM) {
        for (d in 1:nD) {
            alpha[m, d, ] <- as.numeric(
                fit$draws(paste0("alpha_all[", m, ",", d, ",", 1:nT, "]"))
            )
            y[m, d, ] <- as.numeric(
                fit$draws(paste0("y[", m, ",", d, ",", 1:nT, "]"))
            )
        }
    }

    return(list(y = y, alpha = alpha))
}

# ----------------------- Basis Functions---------------------------------

fourier_basis <- function(x,
                          nK,
                          period,
                          intercept = TRUE,
                          derative = 0) {
    # Generates a fourier basis with nK terms with a period for the values in x.

    if (nK %% 2 != 0) {
        stop("nK must be even.")
    }
    omega <- (2 * pi) / period

    if (intercept) {
        basis <- matrix(0, nrow = length(x), ncol = nK + 1)

        ifelse(derative == 0,
            basis[, 1] <- 1,
            basis[, 1] <- 0
        )
        for (k in 1:(nK / 2)) {
            basis[, 2 * k] <- (k * omega)^derative *
                sin((k * omega * x) + (0.5 * derative * pi))
            basis[, 2 * k + 1] <- (k * omega)^derative *
                cos((k * omega * x) + (0.5 * derative * pi))
        }
    } else {
        basis <- matrix(0, nrow = length(x), ncol = nK)
        for (k in 1:(nK / 2)) {
            basis[, 2 * k - 1] <- (k * omega)^derative *
                sin((k * omega * x) + (0.5 * derative * pi))
            basis[, 2 * k] <- (k * omega)^derative *
                cos((k * omega * x) + (0.5 * derative * pi))
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

# ----------------------- Initial Value Functions for Stan Model ---------------
gen_posture <- function(x,
                        Lm,
                        r,
                        theta,
                        degrees = FALSE) {
    # Generate posture given x, Lm, r, theta
    if (degrees) {
        theta <- theta * pi / 180
    }

    nlinks <- 1

    # rotation matrix
    Gam <- vector("list", nlinks)
    for (i in 1:nlinks) {
        Gam[[i]] <- matrix(c(
            cos(theta[i]), sin(theta[i]),
            -sin(theta[i]), cos(theta[i])
        ), 2, 2)
    }

    # compute joint position
    J <- vector("list", nlinks + 1)
    J[[1]] <- matrix(r, nD, 1)
    J[[2]] <- J[[1]] + Gam[[1]] %*% c(0, -Lm)

    # compute alpha
    alpha <- vector(mode = "list", length = nlinks)
    for (i in 1:nlinks) {
        nmarkers <- ncol(x)
        alpha[[i]] <- matrix(NA, 2, nmarkers)
        for (j in 1:nmarkers) {
            x_ij <- x[, j]
            alpha[[i]][, j] <- J[[i]] + Gam[[i]] %*% x_ij
        }
    }
    return(list(r = r, J = J, alpha = alpha))
}

cost <- function(params, y, x, Lm) {
    # Computes the least squares cost of a given pose given observations y
    # and parameters (params)
    # params = vector of parameters (rx, ry, theta1 )
    nlinks <- 1
    r.hat <- params[1:2]
    theta.hat <- params[3:(nlinks + 2)]

    obs.hat <- gen_posture(x, Lm, r = r.hat, theta = theta.hat)
    y.hat <- obs.hat$alpha

    cost <- 0
    nmarkers <- dim(y)[1]
    for (j in 1:nmarkers) {
        diff <- y.hat[[1]][, j] - y[j, ]
        cost <- cost + (sqrt(sum(diff * diff)))
    }
    return(cost)
}

gen_LS_inits <- function(y, x, Lm, inits = c(0, 0, 0)) {
    # Generates initial values via the time-indpenedent LS solution.
    r0 <- inits[1:2]
    theta0 <- inits[3]

    r_inits <<- matrix(NA, nrow = nT, ncol = 2)
    theta_inits <<- rep(NA, nT)
    sigma_inits <<- rep(NA, nT)
    # print(sprintf("Computing time: %.0f", 1))

    # Compute first frame from inits
    temp <- optim(
        par = c(r0, theta0),
        fn = cost,
        y = y[, , 1],
        x = x,
        Lm = lm_true,
        method = "BFGS",
        control = list(trace = FALSE, maxit = 10000, reltol = 1e-8)
    )

    nmarkers <- 0
    obs.hat <- gen_posture(
        x = x, Lm = lm_true,
        r = temp$par[1:2], theta = temp$par[3]
    )
    residual_sum <- 0
    for (j in 1:nM) {
        residual_sum <- residual_sum +
            (y[j, 1, 1] - obs.hat$alpha[[1]][1, j])^2 + # xresidual
            (y[j, 1, 2] - obs.hat$alpha[[1]][2, j])^2 # yresidual
        nmarkers <- nmarkers + 1
    }

    p <- length(temp)
    sigma_inits[1] <- sqrt(residual_sum / (2 * nmarkers))
    r_inits[1, ] <- temp$par[1:2]
    theta_inits[1] <- mod_pi(temp$par[3])


    # compute remaining frames iteratively using previous solution as new initial values
    for (t in 2:nT) {
        # print(sprintf("Computing LS inits at time: %.0f", t))
        temp <- optim(
            par = c(
                r = r_inits[t - 1, ], theta = theta_inits[t - 1]
            ),
            fn = cost,
            y = y[, , t],
            x = x,
            Lm = lm_true,
            method = "BFGS",
            control = list(trace = FALSE, maxit = 10000, reltol = 1e-8)
        )
        r_inits[t, ] <- temp$par[1:2]
        theta_inits[t] <- mod_pi(temp$par[3])

        nmarkers <- 0
        obs.hat <- gen_posture(
            x = x, Lm = lm_true,
            r = temp$par[1:2], theta = temp$par[3]
        )
        residual_sum <- 0
        for (j in 1:nM) {
            residual_sum <- residual_sum +
                (y[j, 1, t] - obs.hat$alpha[[1]][1, j])^2 + # xresidual
                (y[j, 2, t] - obs.hat$alpha[[1]][2, j])^2 # yresidual
            nmarkers <- nmarkers + 1
        }
        p <- length(temp)
        sigma_inits[t] <- sqrt(residual_sum / (2 * nmarkers))
    }
    return(list(
        r_inits = r_inits,
        theta_inits = mod_pi(theta_inits),
        sigma_inits = sigma_inits
    ))
}

pendulum <- function(ts, state, parms) {
    # Pendulum ode for mod 3 initial values
    dydt <- rep(NA, 2)
    dydt[1] <- state[2]
    dydt[2] <- -parms[1] / parms[2] * sin(state[1])
    return(list(dydt))
}

pend_cost <- function(par, y, x) {
    # A LS fit soln to of mod3 used to initialise 
    # par = vector of parameters (rx, ry, theta0, dtheta0, lm)
    ode_soln <- ode(y = par[3:4], times = ts, func = pendulum, parms = c(grav, par[5]))

    err <- rep(NA, nT)
    for (i in 1:nT) {
        obs.hat <- gen_posture(x,
            par[5],
            par[1:2],
            ode_soln[i, 2],
            degrees = FALSE
        )$alpha
        for (j in 1:nM) {
            err[i] <- sum((y[j, , i] - obs.hat[[1]][, j])^2)
        }
    }
    return(mean(err))
}

mod1_inits <- function() {
    # Initialise mod1 with LS solution
    return(list(
        theta = ls_inits$theta_inits,
        lsigma = log(ls_inits$sigma_inits)
    ))
}

mod2_inits <- function() {
    # Initialise mod2 with LS solution

    beta_theta_fit <- lm(ls_inits$theta_inits ~ fbasis - 1)
    beta_lsigma_fit <- lm(log(ls_inits$sigma_inits) ~ fbasis - 1)
    return(list(
        mu_theta = beta_theta_fit$coef[1],
        beta_theta = beta_theta_fit$coef[2:(nK + 1)],
        beta_lsigma = beta_lsigma_fit$coef,
        epsilon_lsigma = sqrt(mean((predict(beta_lsigma_fit) -
            log(ls_inits$sigma_inits))^2)),
        delta_lsigma = rnorm(nT, 0, 0.1)
    ))
}


mod3_inits <- function() {
    # Initialise mod3 with LS solution

    temp <- optim(
        par = c(
            0, 0, -pi / 4, 0, 1
        ),
        fn = pend_cost,
        y = y,
        x = x,
        method = "BFGS",
        control = list(trace = FALSE, maxit = 10000, reltol = 1e-8)
    )

    ode_soln <- ode(y = temp$par[3:4], times = ts, func = pendulum, parms = c(grav, temp$par[5]))

    err <- rep(NA, nT)
    for (i in 1:nT) {
        obs.hat <- gen_posture(x,
            temp$par[5],
            temp$par[1:2],
            ode_soln[i, 2],
            degrees = FALSE
        )$alpha
        for (j in 1:nM) {
            err[i] <- sum((y[j, , i] - obs.hat[[1]][, j])^2)
        }
    }

    return(list(
        theta0 = temp$par[3],
        dtheta0 = temp$par[4],
        lm = temp$par[5],
        lsigma = log(sum(err) / 2 * nT)
    ))
}

# ----------------------- Posterior processing functions -----------------------
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

# ----------------------- Model Performance functions --------- ---------------


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

compute_performance <- function(fit, states_true, sigma_true) {
    # computes performance metrics for a given model fit
    theta_rmse <- rmse(
        states_true[, 1],
        posterior_mean(fit$draws(paste0("theta[", 1:nT, "]")))
    )
    dtheta_rmse <- rmse(
        states_true[, 2],
        posterior_mean(fit$draws(paste0("dtheta[", 1:nT, "]")))
    )
    sigma_rmse <- rmse(
        sigma_true,
        posterior_mean(fit$draws(paste0("sigma[", 1:nT, "]")))
    )

    rmse_tab <- rbind(
        c(round(theta_rmse$rmse, 3), round(theta_rmse$minRMSE, 3), round(theta_rmse$maxRMSE, 3)),
        c(round(dtheta_rmse$rmse, 3), round(dtheta_rmse$minRMSE, 3), round(dtheta_rmse$maxRMSE, 3)),
        c(round(sigma_rmse$rmse * 1000, 3), round(sigma_rmse$minRMSE * 1000, 3), round(sigma_rmse$maxRMSE * 1000, 3))
    )
    rownames(rmse_tab) <- c("theta[rad]", "dtheta[rad/s]", "sigma[mm]")
    colnames(rmse_tab) <- c("mean", "min", "max")

    # CI Width
    theta_CI_Width <- posterior_CI(fit$draws(paste0("theta[", 1:nT, "]")))[, 2] - posterior_CI(fit$draws(paste0("theta[", 1:nT, "]")))[, 1]
    dtheta_CI_Width <- posterior_CI(fit$draws(paste0("dtheta[", 1:nT, "]")))[, 2] - posterior_CI(fit$draws(paste0("dtheta[", 1:nT, "]")))[, 1]
    sigma_CI_Width <- posterior_CI(fit$draws("sigma"))[, 2] - posterior_CI(fit$draws("sigma"))[, 1]

    CI_tab <- round(cbind(
        c(
            mean(theta_CI_Width), mean(dtheta_CI_Width),
            mean(sigma_CI_Width) * 1000
        ),
        c(
            min(theta_CI_Width), min(dtheta_CI_Width),
            min(sigma_CI_Width) * 1000
        ),
        c(
            max(theta_CI_Width), max(dtheta_CI_Width),
            max(sigma_CI_Width) * 1000
        )
    ), 3)

    # Coverage
    theta_coverage <- mean((states_true[, 1] > posterior_CI(fit$draws(paste0("theta[", 1:nT, "]")))[, 1]) & (states_true[, 1] < posterior_CI(fit$draws(paste0("theta[", 1:nT, "]")))[, 2]))
    dtheta_coverage <- mean((states_true[, 2] > posterior_CI(fit$draws(paste0("dtheta[", 1:nT, "]")))[, 1]) & (states_true[, 2] < posterior_CI(fit$draws(paste0("dtheta[", 1:nT, "]")))[, 2]))
    sigma_coverage <- mean((sigma_true > posterior_CI(fit$draws("sigma"))[, 1]) & (sigma_true < posterior_CI(fit$draws("sigma"))[, 2]))
    CI_tab <- cbind(CI_tab, round(c(theta_coverage, dtheta_coverage, sigma_coverage) * 100, 1))

    rownames(CI_tab) <- c("theta[rad]", "dtheta[rad/s]", "sigma[mm]")
    colnames(CI_tab) <- c("mean", "min", "max", "coverage")

    # Sampling efficiency
    theta_sum <- fit$summary(paste0("theta[", 1:nT, "]"))
    dtheta_sum <- fit$summary(paste0("dtheta[", 1:nT, "]"))
    sigma_sum <- fit$summary("sigma")

    ESS_tab <- round(rbind(
        c(mean(theta_sum$ess_bulk), min(theta_sum$ess_bulk), max(theta_sum$ess_bulk)),
        c(mean(dtheta_sum$ess_bulk), min(dtheta_sum$ess_bulk), max(dtheta_sum$ess_bulk)),
        c(mean(sigma_sum$ess_bulk), min(sigma_sum$ess_bulk), max(sigma_sum$ess_bulk))
    ), 0)
    rownames(ESS_tab) <- c("theta[n]", "dtheta[n]", "sigma[n]")
    colnames(ESS_tab) <- c("mean", "min", "max")
    ESS_sec_tab <- ESS_tab / (fit$time()$total)
    rownames(ESS_sec_tab) <- c("theta[n/s]", "dtheta[n/s]", "sigma[n/s]")
    colnames(ESS_sec_tab) <- c("mean", "min", "max")


    Rhat_tab <- round(rbind(
        c(mean(theta_sum$rhat), min(theta_sum$rhat), max(theta_sum$rhat)),
        c(mean(dtheta_sum$rhat), min(dtheta_sum$rhat), max(dtheta_sum$rhat)),
        c(mean(sigma_sum$rhat), min(sigma_sum$rhat), max(sigma_sum$rhat))
    ), 2)
    rownames(Rhat_tab) <- c("theta", "dtheta", "sigma")
    colnames(Rhat_tab) <- c("mean", "min", "max")

    return(list(
        rmse = rmse_tab,
        interval = CI_tab, ess = ESS_tab,
        ess_sec = ESS_sec_tab, rhat = Rhat_tab,
        total_time = fit$time()$total
    ))
}

# ---------------------- Plotting Ultilities -------------------------
# Color plaette for plotting
myColors <<- list(
    grey = rgb(97 / 255, 97 / 255, 97 / 255, 0.1),
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
    r1 = TeX("$r_1 \\,  (m)$"), r2 = TeX("$r_2  \\, (m)$"),
    theta1 = TeX("$\\theta_1 \\, (rad)$"),
    theta2 = TeX("$\\theta_2 \\, (rad)$"),
    theta3 = TeX("$\\theta_3\\, (rad) "),
    sigma = TeX("$\\sigma \\, (mm)$"),
    theta = TeX("$\\theta \\, (rad)$"),
    dtheta = TeX("$\\dot{\\theta} \\, (rad\\cdot s^{-1})$"),
    dthetaBias = TeX("Bias $\\dot{\\theta} \\, (rad\\cdot s^{-1})$"),
    dthetaIntervalWidth = TeX("Credible interval width $\\dot{\\theta} \\, (rad\\cdot s^{-1})$")
)