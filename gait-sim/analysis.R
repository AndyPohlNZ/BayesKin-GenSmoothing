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

# Read true values
true_vals <- readRDS("./true_vals.rdata")
list2env(true_vals, .GlobalEnv)
nT <<- nrow(r_true)
nD <<- ncol(r_true)
nL <<- ncol(theta_true)

# Source statements
source("./library.R")

# Get number of sims to process
nSims <- 100
# --------------------- Process Simulations --------------------------------
mdls <- c("mod1", "mod2")

# preallocate objects to store results
mod1_perf <- vector(mode = "list", length = length(mdls))
mod2_perf <- vector(mode = "list", length = nSims)
thetaPostMeans <- thetaBias <- thetaCoverage <- thetaIntervalWidth <- array(NA, dim = c(length(mdls), nSims, nL, nT))
JAPostMeans <- JABias <- JACoverage <- JAIntervalWidth <- array(NA, dim = c(length(mdls), nSims, 2, nT))
thetaPostCI <- array(NA, dim = c(length(mdls), nSims, 2, nL, nT))
JAPostCI <- array(NA, dim = c(length(mdls), nSims, 2, 2, nT))
rPostMeans <- rBias <- rCoverage <- rIntervalWidth <- array(NA, dim = c(length(mdls), nSims, nD, nT))
rPostCI <- array(NA, dim = c(length(mdls), nSims, 2, nD, nT))

# Process each simulation
for (i in 1:nSims) {
    print(sprintf("Processing simulation %.0f of %.0f", i, nSims))
    for (m in 1:length(mdls)) {
        print(sprintf("---- Processing model: %s", mdls[m]))
        fit <- readRDS(sprintf("./model-fits/%s_%.0f.rds", mdls[m], i))

        # Process performance statistics
        if (m == 1) {
            mod1_perf[[i]] <- compute_performance(r_true, theta_true, sigma_true, fit)
        } else if (m == 2) {
            mod2_perf[[i]] <- compute_performance(r_true, theta_true, sigma_true, fit)
        }

        # Process theta estimates...
        for (l in 1:nL) {
            # theta and estimates
            thetaPostMeans[m, i, l, ] <- posterior_mean(fit$draws(c(paste0("theta[", 1:nT, ",", l, "]"))))
            # theta and dtheta intervals
            thetaPostCI[m, i, 1, l, 1:nT] <- posterior_CI(fit$draws(c(paste0("theta[", 1:nT, ",", l, "]"))))[, 1]
            thetaPostCI[m, i, 2, l, 1:nT] <- posterior_CI(fit$draws(c(paste0("theta[", 1:nT, ",", l, "]"))))[, 2]

            # Bias and coverage
            # theta
            thetaBias[m, i, l, ] <- theta_true[, l] - thetaPostMeans[m, i, l, ]
            thetaCoverage[m, i, l, ] <- (theta_true[, l] > thetaPostCI[m, i, 1, l, ]) & (theta_true[, l] < thetaPostCI[m, i, 2, l, ])
            thetaIntervalWidth[m, i, l, ] <- thetaPostCI[m, i, 2, l, ] - thetaPostCI[m, i, 1, l, ]
        }

        # Process joint angle estimates...
        theta1_draws <- fit$draws(c(paste0("theta[", 1:nT, ",1]")))
        theta2_draws <- fit$draws(c(paste0("theta[", 1:nT, ",2]")))
        theta3_draws <- fit$draws(c(paste0("theta[", 1:nT, ",3]")))
        simsdim <- dim(theta1_draws)
        kja <- aja <- array(NA, dim = simsdim)

        for (ii in 1:simsdim[1]) {
            for (jj in 1:simsdim[2]) {
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
        JAPostMeans[m, i, 1, ] <- posterior_mean(kja)
        JAPostMeans[m, i, 2, ] <- posterior_mean(aja)
        # theta and dtheta intervals
        JAPostCI[m, i, 1, 1, 1:nT] <- posterior_CI(kja)[, 1]
        JAPostCI[m, i, 2, 1, 1:nT] <- posterior_CI(kja)[, 2]
        JAPostCI[m, i, 1, 2, 1:nT] <- posterior_CI(aja)[, 1]
        JAPostCI[m, i, 2, 2, 1:nT] <- posterior_CI(aja)[, 2]
        # Bias and coverage
        kja_true <- knee_joint_angle(theta_true[, 1], theta_true[, 2])
        aja_true <- ankle_joint_angle(theta_true[, 2], theta_true[, 3])

        JABias[m, i, 1, ] <- kja_true - JAPostMeans[m, i, 1, ]
        JABias[m, i, 2, ] <- aja_true - JAPostMeans[m, i, 2, ]
        JACoverage[m, i, 1, ] <- (kja_true > JAPostCI[m, i, 1, 1, ]) & (kja_true < JAPostCI[m, i, 2, 1, ])
        JACoverage[m, i, 2, ] <- (aja_true > JAPostCI[m, i, 1, 2, ]) & (aja_true < JAPostCI[m, i, 2, 2, ])
        JAIntervalWidth[m, i, 1, ] <- JAPostCI[m, i, 2, 1, ] - JAPostCI[m, i, 1, 1, ]
        JAIntervalWidth[m, i, 2, ] <- JAPostCI[m, i, 2, 2, ] - JAPostCI[m, i, 1, 2, ]


        # process r estimates
        for (d in 1:nD) {
            # theta and estimates
            rPostMeans[m, i, d, ] <- posterior_mean(fit$draws(c(paste0("r[", 1:nT, ",", d, "]"))))
            # theta and dtheta intervals
            rPostCI[m, i, 1, d, 1:nT] <- posterior_CI(fit$draws(c(paste0("r[", 1:nT, ",", d, "]"))))[, 1]
            rPostCI[m, i, 2, d, 1:nT] <- posterior_CI(fit$draws(c(paste0("r[", 1:nT, ",", d, "]"))))[, 2]

            # Bias and coverage
            # theta
            rBias[m, i, d, ] <- r_true[, d] - rPostMeans[m, i, d, ]
            rCoverage[m, i, d, ] <- (r_true[, d] > rPostCI[m, i, 1, d, ]) & (r_true[, d] < rPostCI[m, i, 2, d, ])
            rIntervalWidth[m, i, d, ] <- rPostCI[m, i, 2, d, ] - rPostCI[m, i, 1, d, ]
        }
    }
}

# Save results
save.image(file = "./model-fits/processedFiles.rda")

# --------------------- Examinine sampling time --------------------------------
mean(sapply(mod2_perf, function(x) {
    return(x$total_time)
}))

# --------------------- Generate plots referenced in the manuscript-------------
load("./model-fits/processedFiles.rda")
nSims <- 100
model_titles <- c("GAIT_MOD1: Time Independent BIK", "GAIT_MOD2: Generalised Smoothing BIK")
modelColors <- list(
    rgb(247 / 255, 129 / 255, 191 / 255, 0.2),
    rgb(55 / 255, 78 / 255, 163 / 255, 0.2),
    rgb(77 / 255, 175 / 255, 74 / 255, 0.2)
)


# Raw estimates - Joint angles.
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Gait_JA.pdf",
    paper = "letter"
)
kja_true <- knee_joint_angle(theta_true[, 1], theta_true[, 2])
aja_true <- ankle_joint_angle(theta_true[, 2], theta_true[, 3])
ja_true <- rbind(kja_true, aja_true)
par(mfrow = c(2, 2))
for (a in 1:2) { # TODO change once aja function fixed
    for (i in 1:length(mdls)) {
        plot(NA,
            xlim = range(ts), ylim = c(-0.5, pi / 2),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$JA[a], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(a == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(a == 1, model_titles[i], ""), line = -1, cex.main = 1.5)

        for (j in sample(1:nSims, size = 20)) {
            lines(ts, mod_pi(JAPostMeans[i, j, a, ]), col = myColors$grey)
        }
        lines(ts, ja_true[a, ], col = myColors$orange, lty = 1, lwd = 1.5)
    }
}
dev.off()

# JA Bias
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Bias_JA.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (a in 1:2) {
    for (i in 1:length(mdls)) {
        # bias
        plot(NA,
            xlim = range(ts), ylim = c(-0.1, 0.1),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$JA_bias[a], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(a == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(a == 1, model_titles[i], ""), line = -1, cex.main = 1.5)


        segments(
            x0 = ts, y0 = apply(JABias[i, , a, ], 2, function(x) {
                quantile(x, probs = c(0.1))
            }),
            x1 = ts, y1 = apply(JABias[i, , a, ], 2, function(x) {
                quantile(x, probs = c(0.9))
            }),
            col = "gray77"
        )
        points(ts, apply(JABias[i, , a, ], 2, mean), pch = 16, cex = 0.3)
        lines(ts, rep(0, nT), col = myColors$orange)
    }
}
dev.off()

# RMSE

pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/RMSE_JA.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (a in 1:2) {
    for (i in 1:length(mdls)) {
        rmse <- sqrt(colMeans(JABias[i, , a, ]^2))
        print(sprintf("Model %.0f, angle %.0f: RMSE = %.3f", i, a, mean(rmse)))
        plot(NA,
            xlim = range(ts), ylim = c(0, 0.1),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$JA_RMSE[a], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(a == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(a == 1, model_titles[i], ""), line = -1, cex.main = 1.5)

        lines(ts, rmse, col = "gray77")
    }
}
dev.off()



# JA Coverage
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Coverage_JA.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (a in 1:2) {
    for (i in 1:length(mdls)) {
        # Coverage
        plot(NA,
            xlim = range(ts), ylim = c(0.5, 1),
            main = "", ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$JA_coverage[a], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(a == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(a == 1, model_titles[i], ""), line = -1, cex.main = 1.5)

        points(ts, apply(JACoverage[i, , a, ], 2, mean), pch = 16, cex = 0.5)
        lines(ts, rep(0.95, nT), col = myColors$orange)
    }
}
dev.off()

# JA Interval width
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/CIWidth_JA.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
ylims <- list(c(0, 0.2), c(0, 0.5))
for (a in 1:2) {
    for (i in 1:length(mdls)) {
        plot(NA,
            xlim = range(ts), ylim = ylims[[a]],
            main = "", ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$JA_intervalWidth[a], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(a == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(a == 1, model_titles[i], ""), line = -1, cex.main = 1.5)
        segments(
            x0 = ts, y0 = apply(JAIntervalWidth[i, , a, ], 2, function(x) {
                quantile(x, probs = c(0.1))
            }),
            x1 = ts, y1 = apply(JAIntervalWidth[i, , a, ], 2, function(x) {
                quantile(x, probs = c(0.9))
            }),
            col = "gray77"
        )
        points(ts, apply(JAIntervalWidth[i, , a, ], 2, mean), pch = 16, cex = 0.3)
    }
}
dev.off()


mean(apply(JAIntervalWidth[1, , 2, ], 2, mean))
mean(apply(JAIntervalWidth[2, , 2, ], 2, mean))
min(sapply(mod1_perf, function(x) {
    return(min(x$ess[, 1]))
}))

# ----------------------- Model Parameters -----------------------


# Raw estimates - theta
pdf(
    width = 7 * 1, height = 7 * 3 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Gait_theta.pdf",
    paper = "letter"
)
par(mfrow = c(3, 2))
ylims <- list(c(-2.2, -1), c(-3.2, -1.2), c(-2, 0))
for (l in 1:nL) {
    for (i in 1:length(mdls)) {
        plot(NA,
            xlim = range(ts), ylim = ylims[[l]],
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )

        title(ylab = parmLabs$theta[l], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(l == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(l == 1, model_titles[i], ""), line = -1, cex.main = 1.5)

        for (j in 1:(nSims / 2)) {
            lines(ts, thetaPostMeans[i, j, l, ], col = myColors$grey)
        }
        lines(ts, theta_true[, l], col = myColors$orange, lty = 1, lwd = 1.5)
    }
}
dev.off()


# Raw estiamtes - r
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Gait_r.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (d in 1:nD) {
    for (i in 1:length(mdls)) {
        plot(NA,
            xlim = range(ts), ylim = range(rPostMeans[i, , d, ]),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$r[d], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(d == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(d == 1, model_titles[i], ""), line = -1, cex.main = 1.5)
        for (j in 1:(nSims / 2)) {
            lines(ts, rPostMeans[i, j, d, ], col = myColors$grey)
        }
        lines(ts, r_true[, d], col = myColors$orange, lty = 1, lwd = 1.5)
    }
}
dev.off()


# Bias - theta
pdf(
    width = 7 * 1, height = 3 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Bias_theta.pdf",
    paper = "letter"
)
par(mfrow = c(3, 2))
for (l in 1:nL) {
    for (i in 1:length(mdls)) {
        # bias
        plot(NA,
            xlim = range(ts), ylim = c(-0.1, 0.1),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$biasTheta[l], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(l == 3, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(l == 1, model_titles[i], ""), line = -1, cex.main = 1.5)
        segments(
            x0 = ts, y0 = apply(thetaBias[i, , l, ], 2, function(x) {
                quantile(x, probs = c(0.1))
            }),
            x1 = ts, y1 = apply(thetaBias[i, , l, ], 2, function(x) {
                quantile(x, probs = c(0.9))
            }),
            col = "gray77"
        )
        points(ts, apply(thetaBias[i, , l, ], 2, mean), pch = 16, cex = 0.3)
        lines(ts, rep(0, nT), col = myColors$orange)
    }
}
dev.off()

# Bias - r
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Bias_r.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (d in 1:nD) {
    for (i in 1:length(mdls)) {
        # bias
        plot(NA,
            xlim = range(ts), ylim = c(-0.03, 0.03),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$biasR[d], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(d == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(d == 1, model_titles[i], ""), line = -1, cex.main = 1.5)
        segments(
            x0 = ts, y0 = apply(rBias[i, , d, ], 2, function(x) {
                quantile(x, probs = c(0.1))
            }),
            x1 = ts, y1 = apply(rBias[i, , d, ], 2, function(x) {
                quantile(x, probs = c(0.9))
            }),
            col = "gray77"
        )

        points(ts, apply(rBias[i, , d, ], 2, mean), pch = 16, cex = 0.3)
        lines(ts, rep(0, nT), col = myColors$orange)
    }
}
dev.off()

# RMSE theta

pdf(
    width = 7 * 1, height = 3 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/RMSE_theta.pdf",
    paper = "letter"
)
par(mfrow = c(3, 2))
for (l in 1:3) {
    for (i in 1:length(mdls)) {
        rmse <- sqrt(colMeans(thetaBias[i, , l, ]^2))
        print(sprintf("Model %.0f, angle %.0f: RMSE = %.3f", i, l, mean(rmse)))
        plot(NA,
            xlim = range(ts), ylim = c(0, 0.1),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$RMSETheta[l], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(l == 3, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(l == 1, model_titles[i], ""), line = -1, cex.main = 1.5)

        lines(ts, rmse, col = "gray77")
    }
}
dev.off()

pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/RMSE_r.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (d in 1:2) {
    for (i in 1:length(mdls)) {
        rmse <- sqrt(colMeans(rBias[i, , d, ]^2))
        print(sprintf("Model %.0f, dimension %.0f: RMSE = %.3f", i, d, mean(rmse)))
        plot(NA,
            xlim = range(ts), ylim = c(0, 0.015),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$RMSER[d], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(d == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(d == 1, model_titles[i], ""), line = -1, cex.main = 1.5)

        lines(ts, rmse, col = "gray77")
    }
}
dev.off()




# Coverage - theta
pdf(
    width = 7 * 1, height = 3 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Coverage_theta.pdf",
    paper = "letter"
)
par(mfrow = c(3, 2))
for (l in 1:nL) {
    for (i in 1:length(mdls)) {
        # Coverage
        plot(NA,
            xlim = range(ts), ylim = c(0.5, 1),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$coverageTheta[l], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(l == 3, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(l == 1, model_titles[i], ""), line = 1, cex.main = 1.5)
        points(ts, apply(thetaCoverage[i, , l, ], 2, mean), pch = 16, cex = 0.5)
        lines(ts, rep(0.95, nT), col = myColors$orange)
    }
}
dev.off()

# Coverage - r
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/Coverage_r.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (d in 1:nD) {
    for (i in 1:length(mdls)) {
        # Coverage
        plot(NA,
            xlim = range(ts), ylim = c(0.5, 1),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$coverageR[d], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(d == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(d == 1, model_titles[i], ""), line = 1, cex.main = 1.5)
        points(ts, apply(rCoverage[i, , d, ], 2, mean), pch = 16, cex = 0.5)
        lines(ts, rep(0.95, nT), col = myColors$orange)
    }
}
dev.off()


# Interval Width - theta
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/CIWidth_theta.pdf",
    paper = "letter"
)
par(mfrow = c(3, 2))
for (l in 1:nL) {
    for (i in 1:length(mdls)) {
        plot(NA,
            xlim = range(ts), ylim = c(0, 0.35),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$intervalWidthTheta[l], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(l == 3, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(l == 1, model_titles[i], ""), line = -1, cex.main = 1.5)
        segments(
            x0 = ts, y0 = apply(thetaIntervalWidth[i, , l, ], 2, function(x) {
                quantile(x, probs = c(0.1))
            }),
            x1 = ts, y1 = apply(thetaIntervalWidth[i, , l, ], 2, function(x) {
                quantile(x, probs = c(0.9))
            }),
            col = "gray77"
        )
        points(ts, apply(thetaIntervalWidth[i, , l, ], 2, mean), pch = 16, cex = 0.3)
    }
}
dev.off()

# Interval Width - r
pdf(
    width = 7 * 1, height = 2 * 7 * 0.4, pointsize = 12 * 0.5,
    file = "./figures/CIWidth_r.pdf",
    paper = "letter"
)
par(mfrow = c(2, 2))
for (d in 1:nD) {
    for (i in 1:length(mdls)) {
        plot(NA,
            xlim = range(ts), ylim = c(0, 0.06),
            main = "",
            ylab = "", xlab = "",
            bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5
        )
        title(ylab = parmLabs$intervalWidthR[d], line = 2, cex.lab = 1.5)
        title(xlab = ifelse(d == 2, "Time (s)", ""), line = 2.2, cex.lab = 1.5)
        title(main = ifelse(d == 1, model_titles[i], ""), line = -1, cex.main = 1.5)
        segments(
            x0 = ts, y0 = apply(rIntervalWidth[i, , d, ], 2, function(x) {
                quantile(x, probs = c(0.1))
            }),
            x1 = ts, y1 = apply(rIntervalWidth[i, , d, ], 2, function(x) {
                quantile(x, probs = c(0.9))
            }),
            col = "gray77"
        )
        points(ts, apply(rIntervalWidth[i, , d, ], 2, mean), pch = 16, cex = 0.3)
    }
}
dev.off()