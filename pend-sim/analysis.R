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


# --------------------- Perliminaries -------------------------------------

# Read true values
true_vals <- readRDS("./true_vals.rds")
list2env(true_vals, .GlobalEnv)

# Source statements
source("./library.R")

# Get number of sims to process
nSims <- 100
# --------------------- Process Simulations --------------------------------

mdls <- c("mod1", "mod2", "mod3")
mod_perf <- vector(mode = "list", length = length(mdls))
mod1_perf <- vector(mode = "list", length = nSims)
mod2_perf <- list()
mod3_perf <- list()
thetaPostMeans <- thetaBias <- thetaCoverage <- thetaIntervalWidth <- array(NA, dim = c(length(mdls), nSims, nT))
dthetaPostMeans <- dthetaBias <- dthetaCoverage <- dthetaIntervalWidth <- array(NA, dim = c(length(mdls), nSims, nT))
thetaPostCI <- array(NA, dim = c(length(mdls), nSims, 2, nT))
dthetaPostCI <- array(NA, dim = c(length(mdls), nSims, 2, nT))
for (i in 1:nSims) {
    print(sprintf("Processing simulation %.0f of %.0f", i, nSims))
    for (m in 1:length(mdls)) {
        print(sprintf("---- Processing model: %s", mdls[m]))
        fit <- readRDS(sprintf("./model-fits/%s_fit_%.0f.rds", mdls[m], i))

        # Process performance statistics
        if (m == 1) {
            mod1_perf[[i]] <- compute_performance(fit, states_true, sigma_true)
        } else if (m == 2) {
            mod2_perf[[i]] <- compute_performance(fit, states_true, sigma_true)
        } else if (m == 3) {
            mod3_perf[[i]] <- compute_performance(fit, states_true, sigma_true)
        }

        theta_samps <- fit$draws(c("theta"))
        dtheta_samps <- fit$draws(c("dtheta"))

        # theta and dtheta estimates
        thetaPostMeans[m, i, ] <- posterior_mean(theta_samps)
        dthetaPostMeans[m, i, ] <- posterior_mean(dtheta_samps)

        # theta and dtheta intervals
        thetaPostCI[m, i, 1, 1:nT] <- posterior_CI(theta_samps)[, 1]
        thetaPostCI[m, i, 2, 1:nT] <- posterior_CI(theta_samps)[, 2]

        dthetaPostCI[m, i, 1, 1:nT] <- posterior_CI(dtheta_samps)[, 1]
        dthetaPostCI[m, i, 2, 1:nT] <- posterior_CI(dtheta_samps)[, 2]

        # Bias and coverage
        # theta
        thetaBias[m, i, ] <- states_true[, 1] - thetaPostMeans[m, i, ]
        thetaCoverage[m, i, ] <- (states_true[, 1] > thetaPostCI[m, i, 1, ]) & (states_true[, 1] < thetaPostCI[m, i, 2, ])
        thetaIntervalWidth[m, i, ] <- thetaPostCI[m, i, 2, ] - thetaPostCI[m, i, 1, ]
        # dtheta
        dthetaBias[m, i, ] <- states_true[, 2] - dthetaPostMeans[m, i, ]
        dthetaCoverage[m, i, ] <- (states_true[, 2] > dthetaPostCI[m, i, 1, ]) & (states_true[, 2] < dthetaPostCI[m, i, 2, ])
        dthetaIntervalWidth[m, i, ] <- dthetaPostCI[m, i, 2, ] - dthetaPostCI[m, i, 1, ]
    }
}

# Compute RMSE
model_rmse_tmp <- array(NA, dim = c(length(mdls), nSims, nT))
model_rmse <- array(NA, dim = c(length(mdls), nT))
model_rmse_tmp_dtheta <- array(NA, dim = c(length(mdls), nSims, nT))
model_rmse_dtheta <- array(NA, dim = c(length(mdls), nT))
for (m in 1:length(mdls)) {
    for (i in 1:nSims) {
        model_rmse_tmp[m, i, ] <- thetaBias[m, i, ]^2
        model_rmse_tmp_dtheta[m, i, ] <- dthetaBias[m, i, ]^2
    }
    model_rmse[m, ] <- sqrt(apply(model_rmse_tmp[m, , ], 2, mean))
    model_rmse_dtheta[m, ] <- sqrt(apply(model_rmse_tmp_dtheta[m, , ], 2, mean))
}

# Save results
save.image(file = "./model-fits/processedFiles.rda")

# --------------------- Plot Results -------------------------------------------

# If pre processed load results:
# load("./model-fits/processedFiles.rda")

nSims <- 100
model_titles <- c("PEND_MOD1: Time Independent BIK", "PEND_MOD2: Generalised Smoothing BIK", "PEND_MOD3: 'True' ODE model")
modelColors <- list(
    rgb(247 / 255, 129 / 255, 191 / 255, 1),
    rgb(55 / 255, 78 / 255, 163 / 255, 1),
    rgb(77 / 255, 175 / 255, 74 / 255, 1)
)

# --------------------- theta---------------------------------------
# Raw estimates
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_theta.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    plot(NA,
        xlim = range(ts), ylim = range(thetaPostMeans[, , ]),
        main = "",
        ylab = "", xlab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
    )
    title(ylab = parmLabs$theta, line = 2, cex.lab = 2)
    title(main = model_titles[i], line = -1, cex.main = 1.5)

    for (j in sample(1:nSims, size = 20)) {
        lines(ts, thetaPostMeans[i, j, ], col = "grey")
    }
    lines(ts, states_true[, 1], col = myColors$orange, lty = 1, lwd = 0.9)
}
dev.off()

# Bias
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_Bias.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    # bias
    plot(NA,
        xlim = range(ts), ylim = c(-0.08, 0.08),
        ylab = "", xlab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = "Bias (rad)", , line = 2, cex.lab = 2)

    segments(
        x0 = ts, y0 = apply(thetaBias[i, , ], 2, function(x) {
            quantile(x, probs = c(0.1))
        }),
        x1 = ts, y1 = apply(thetaBias[i, , ], 2, function(x) {
            quantile(x, probs = c(0.9))
        }),
        col = "gray77"
    )


    points(ts, apply(thetaBias[i, , ], 2, mean), pch = 16, cex = 0.3)
    lines(ts, rep(0, nT), col = myColors$orange)
}
dev.off()

# RMSE
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_RMSE.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    # bias
    plot(NA,
        xlim = range(ts), ylim = c(0, max(model_rmse)),
        ylab = "", xlab = "Time (s)",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = "RMSE (rad)", , line = 2, cex.lab = 2)

    lines(ts, model_rmse[i, ], col = "grey77")
}
dev.off()

# Coverage
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_Coverage.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    # Coverage
    plot(NA,
        xlim = range(ts), ylim = c(0.5, 1),
        xlab = "", ylab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )

    title(ylab = "Coverage (%)", , line = 2, cex.lab = 2)
    title(main = model_titles[i], line = 0, cex.main = 1.5)

    points(ts, apply(thetaCoverage[i, , ], 2, mean), pch = 16, cex = 0.5)
    lines(ts, rep(0.95, nT), col = myColors$orange)
}
dev.off()

# Interval Width
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_CIWidth.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    plot(NA,
        xlim = range(ts), ylim = c(0, 0.25),
        ylab = "", xlab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = "Credible Interval width (rad)", , line = 2, cex.lab = 2)
    title(xlab = "Time (s)", , line = 2, cex.lab = 2)

    segments(
        x0 = ts, y0 = apply(thetaIntervalWidth[i, , ], 2, function(x) {
            quantile(x, probs = c(0.1))
        }),
        x1 = ts, y1 = apply(thetaIntervalWidth[i, , ], 2, function(x) {
            quantile(x, probs = c(0.9))
        }),
        col = "gray77"
    )
    points(ts, apply(thetaIntervalWidth[i, , ], 2, mean), pch = 16, cex = 0.3)
}
dev.off()

mean(apply(thetaIntervalWidth[3, , ], 2, mean))
# --------------------- dtheta---------------------------------------

# Raw estimates
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_dtheta.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
ylims <- list(c(-10, 10), c(-2, 2), c(-2, 2))

for (i in 1:length(mdls)) {
    plot(NA,
        xlim = range(ts), ylim = ylims[[i]],
        main = "",
        ylab = "", xlab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = parmLabs$dtheta, line = 2, cex.lab = 2)
    title(main = model_titles[i], line = 1, cex.main = 1.5)

    for (j in sample(1:nSims, size = 20)) {
        lines(ts, dthetaPostMeans[i, j, ], col = myColors$grey)
    }
    lines(ts, states_true[, 2], col = myColors$orange, lty = 1, lwd = 1.5)
}
dev.off()


# Bias
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_dtheta_Bias.pdf",
    paper = "letter"
)
ylims <- list(c(-10, 10), c(-1, 1), c(-0.1, 0.1))
par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    # bias
    plot(NA,
        xlim = range(ts), ylim = ylims[[i]],
        ylab = "",
        xlab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = TeX("Bias: $(rad\\cdot s^{-1})$"), , line = 2, cex.lab = 2)

    segments(
        x0 = ts, y0 = apply(dthetaBias[i, , ], 2, function(x) {
            quantile(x, probs = c(0.1))
        }),
        x1 = ts, y1 = apply(dthetaBias[i, , ], 2, function(x) {
            quantile(x, probs = c(0.9))
        }),
        col = "gray77"
    )

    points(ts, apply(dthetaBias[i, , ], 2, mean), pch = 16, cex = 0.3)
    lines(ts, rep(0, nT), col = myColors$orange)
}
dev.off()

# RMSE
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_dtheta_RMSE.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
ylims <- list(c(0, 10), c(0, 0.2), c(0, 0.1))
for (i in 1:length(mdls)) {
    plot(NA,
        xlim = range(ts), ylim = ylims[[i]],
        ylab = "", xlab = "Time (s)",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = "RMSE (rad)", , line = 2, cex.lab = 2)

    lines(ts, model_rmse_dtheta[i, ], col = "grey77")
}
dev.off()

# Coverage
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_dtheta_Coverage.pdf",
    paper = "letter"
)

par(mfrow = c(1, 3))
for (i in 1:length(mdls)) {
    # Coverage
    plot(NA,
        xlim = range(ts), ylim = c(0.5, 1),
        xlab = "", ylab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = "Coverage (%)", , line = 2, cex.lab = 2)
    title(main = model_titles[i], line = 0, cex.main = 1.5)
    points(ts, apply(dthetaCoverage[i, , ], 2, mean), pch = 16, cex = 0.5)
    lines(ts, rep(0.95, nT), col = myColors$orange)
}
dev.off()

# Interval Width
pdf(
    width = 7 * 1, height = 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Pend_dtheta_CIWidth.pdf",
    paper = "letter"
)
par(mfrow = c(1, 3))
ylims <- list(c(0, 30), c(0, 1), c(0, 0.1))
for (i in 1:length(mdls)) {
    plot(NA,
        xlim = range(ts), ylim = ylims[[i]],
        ylab = "",
        xlab = "",
        bty = "n", cex.axis = 1.5, cex.lab = 1.5
    )
    title(ylab = TeX("Credible interval width $(rad\\cdot s^{-1})$"), line = 2, cex.lab = 2)
    title(xlab = "Time (s)", , line = 2, cex.lab = 2)

    segments(
        x0 = ts, y0 = apply(dthetaIntervalWidth[i, , ], 2, function(x) {
            quantile(x, probs = c(0.1))
        }),
        x1 = ts, y1 = apply(dthetaIntervalWidth[i, , ], 2, function(x) {
            quantile(x, probs = c(0.9))
        }),
        col = "gray77"
    )
    points(ts, apply(dthetaIntervalWidth[i, , ], 2, mean), pch = 16, cex = 0.3)
}
dev.off()