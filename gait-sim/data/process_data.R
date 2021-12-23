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

# ------------------------------------------------------------------------------

# Process ns strides from raw kinematic data located in data directory
## Input
# datadir = 'FILENAME'
# ns = number of strides
# first_stide = first stride

## Output
# x - static trial
# y - observed marker locations
# ts - frame times.


process_data <- function(outputdir = "./data",
                         datadir = "./data",
                         nS = 1,
                         first_stride = 5,
                         plot = FALSE) {
    # Basic functions
    gen_rotation_matrix <- function(theta) {
        # Generates a rotation matrix for angle theta specified in radians

        if (theta > 0) {
            return(rbind(c(cos(theta), -sin(theta)), c(sin(theta), cos(theta))))
        }
    }
    gen_vec_angle <- function(v) {
        # Generates the angle a vector v makes with the x axis
        angle <- abs(atan2(v[2], v[1]))
        return(angle)
    }

    # ---------------------- Data Loading/Processing --------------------------------

    # Z-Y = sagittal plane...
    # load data for subject.
    static <- read.table(sprintf("%s/s2_staticTrial.csv", datadir), sep = ",", header = TRUE)
    events <- read.table(sprintf("%s//s2_events.csv", datadir), sep = ",", header = TRUE)
    events <- events[complete.cases(events), ]
    dynamic <- read.table(sprintf("%s//s2_dynamicTrial.csv", datadir), sep = ",", header = TRUE)

    # Re-align coordinate system to sensible x, y z system.
    static$Z <- -static$Z # reflect
    static[, c("X", "Y", "Z")] <- static[, c("X", "Y", "Z")] / 1000 # change to m

    dynamic$Z <- -dynamic$Z # reflect
    dynamic[, c("X", "Y", "Z")] <- dynamic[, c("X", "Y", "Z")] / 1000 # change to m

    # Comptue Joint Centers
    HJC <- static[static$marker == "RHip", c("Z", "Y")]
    KJC <- static[static$marker == "LatKnee", c("Z", "Y")] + (static[static$marker == "MedKnee", c("Z", "Y")] - static[static$marker == "LatKnee", c("Z", "Y")]) / 2
    AJC <- static[static$marker == "LatAnkle", c("Z", "Y")] + (static[static$marker == "MedAnkle", c("Z", "Y")] - static[static$marker == "LatAnkle", c("Z", "Y")]) / 2
    TOE <- static[static$marker == "Toe", c("Z", "Y")]

    # Compute LCSs
    thigh <- rbind(
        static[static$marker == "Thigh1", c("Z", "Y")],
        static[static$marker == "Thigh2", c("Z", "Y")],
        static[static$marker == "Thigh3", c("Z", "Y")],
        static[static$marker == "Thigh4", c("Z", "Y")]
    )


    shank <- rbind(
        static[static$marker == "Shank1", c("Z", "Y")],
        static[static$marker == "Shank2", c("Z", "Y")],
        static[static$marker == "Shank3", c("Z", "Y")],
        static[static$marker == "Shank4", c("Z", "Y")]
    )


    foot <- rbind(
        static[static$marker == "Foot1", c("Z", "Y")],
        static[static$marker == "Foot2", c("Z", "Y")],
        static[static$marker == "Foot3", c("Z", "Y")],
        static[static$marker == "Foot4", c("Z", "Y")]
    )

    if (plot) {
        # plot static trial
        plot(rbind(HJC, KJC),
            type = "o", pch = 16, col = rgb(18 / 255, 78 / 255, 120 / 255, 1),
            xlim = c(-0.6, 0.0), ylim = c(0, 1), asp = 1,
            xlab = "x (m)", ylab = "y (m)",
            bty = "n"
        )
        points(thigh, pch = 4, col = rgb(18 / 255, 78 / 255, 120 / 255, 1))
        lines(rbind(KJC, AJC), type = "o", pch = 16, col = rgb(76 / 255, 159 / 255, 112 / 255, 1))
        points(shank, pch = 4, col = rgb(76 / 255, 159 / 255, 112 / 255, 1))
        lines(rbind(AJC, TOE), type = "o", pch = 16, col = rgb(217 / 255, 3 / 255, 104 / 255, 1))
        points(foot, pch = 4, col = rgb(217 / 255, 3 / 255, 104 / 255, 1))
    }

    # Realign to the horizontal static trial expected by Bayesian Model
    thigh_lcs_origin <- c(0, 0)
    thigh_lcs_termination <- as.numeric(KJC - HJC)
    angle <- gen_vec_angle(thigh_lcs_termination)
    rot.mat <- gen_rotation_matrix(angle)
    thigh_lcs_termination <- rot.mat %*% thigh_lcs_termination
    thigh_lcs <- do.call(rbind, apply(thigh, 1, function(x) {
        return(x - HJC)
    }))
    thigh_lcs <- apply(thigh_lcs, 1, function(x) {
        return(rot.mat %*% x)
    })

    shank_lcs_origin <- c(0, 0)
    shank_lcs_termination <- as.numeric(AJC - KJC)
    angle <- gen_vec_angle(shank_lcs_termination)
    rot.mat <- gen_rotation_matrix(angle)
    shank_lcs_termination <- rot.mat %*% shank_lcs_termination
    shank_lcs <- do.call(rbind, apply(shank, 1, function(x) {
        return(x - KJC)
    }))
    shank_lcs <- apply(shank_lcs, 1, function(x) {
        return(rot.mat %*% x)
    })

    foot_lcs_origin <- c(0, 0)
    foot_lcs_termination <- as.numeric(TOE - AJC)
    angle <- gen_vec_angle(foot_lcs_termination)
    rot.mat <- gen_rotation_matrix(angle)
    foot_lcs_termination <- rot.mat %*% foot_lcs_termination
    foot_lcs <- do.call(rbind, apply(foot, 1, function(x) {
        return(x - AJC)
    }))
    foot_lcs <- apply(foot_lcs, 1, function(x) {
        return(rot.mat %*% x)
    })

    # Add row/colnames
    rownames(thigh_lcs) <- c("X", "Y")
    colnames(thigh_lcs) <- c("thigh1", "thigh2", "thigh3", "thigh4")
    rownames(shank_lcs) <- c("X", "Y")
    colnames(shank_lcs) <- c("shank1", "shank2", "shank3", "shank4")
    rownames(foot_lcs) <- c("X", "Y")
    colnames(foot_lcs) <- c("foot1", "foot2", "foot3", "foot4")

    if (plot) {
        # Plot 'x'
        plot(rbind(thigh_lcs_origin, t(thigh_lcs_termination)),
            type = "o", pch = 16,
            xlim = c(-1, 1), ylim = c(-.8, 1), asp = 1,
            main = "Thigh LCS", xlab = "x", ylab = "y"
        )
        points(t(thigh_lcs), pch = 4)


        plot(rbind(shank_lcs_origin, t(shank_lcs_termination)),
            type = "o", pch = 16,
            xlim = c(-1, 1), ylim = c(-.8, 1), asp = 1,
            main = "Shank LCS", xlab = "x", ylab = "y"
        )
        points(t(shank_lcs), pch = 4)


        plot(rbind(foot_lcs_origin, t(foot_lcs_termination)),
            type = "o", pch = 16,
            xlim = c(-1, 1), ylim = c(-.8, 1), asp = 1,
            main = "Foot LCS", xlab = "x", ylab = "y"
        )
        points(t(foot_lcs), pch = 4)
    }

    # Compute Length of each segment
    thigh_L <- sqrt(sum(thigh_lcs_termination - thigh_lcs_origin)^2)
    shank_L <- sqrt(sum(shank_lcs_termination - shank_lcs_origin)^2)
    foot_L <- sqrt(sum(foot_lcs_termination - foot_lcs_origin)^2)


    # Generalte links format
    link1 <- list(
        length = thigh_L,
        x = thigh_lcs
    )
    link2 <- list(
        length = shank_L,
        x = shank_lcs
    )
    link3 <- list(
        length = foot_L,
        x = foot_lcs
    )
    links <- list(thigh = link1, shank = link2, foot = link3)

    # ----------------------- Prepare constants -------------------------------------------
    # Generate constants
    nL <<- length(links)
    nM <<- max(sapply(links, function(x) {
        ncol(x$x)
    }))
    nD <<- 2
    dT <<- 0.005 # sampling period
    sampfreq <<- 1 / dT # sampling frequency

    # Trim dynamic and times to nS complete gait cycles (FS to FS) offset by 'first_stride'
    events <- events * dT # convert events from index to time
    mint <- events[first_stride, 1] # starting with the 'first_stride'
    maxt <- events[nS + first_stride, 1]
    dynamic <- dynamic[(dynamic$time >= mint & dynamic$time <= maxt), ]
    dynamic$time <- dynamic$time - min(dynamic$time) # reset start time at 0 sec.
    events <- events[first_stride:(nS + first_stride), ] - mint

    maxtime <<- max(dynamic$time)
    nT <<- maxtime / dT + 1 # number of time points


    #### Static trial
    x <- array(data = NA, dim = c(nL, nM, nD)) # note dimension so that x[1,,1] is a vector...
    xl <- lapply(links, function(x) {
        x$x
    })
    x[1, 1:nM, 1:nD] <- t(xl[[1]])
    x[2, 1:nM, 1:nD] <- t(xl[[2]])
    x[3, 1:nM, 1:nD] <- t(xl[[3]])

    # ----------------------- Prepare observations -------------------------------------------
    ts <- seq(0, maxtime, by = dT)
    y <- array(NA, dim = c(nL, max(nM), nT, nD))
    for (t in 1:nT) {
        tidx <- which(abs(dynamic$time - ts[t]) < 1e-12)
        dynamic_t <- dynamic[tidx, ]

        y[1, , t, ] <- as.matrix(rbind(
            dynamic_t[dynamic_t$marker == "Thigh1", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Thigh2", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Thigh3", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Thigh4", c("Z", "Y")]
        ))

        y[2, , t, ] <- as.matrix(rbind(
            dynamic_t[dynamic_t$marker == "Shank1", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Shank2", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Shank3", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Shank4", c("Z", "Y")]
        ))

        y[3, , t, ] <- as.matrix(rbind(
            dynamic_t[dynamic_t$marker == "Foot1", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Foot2", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Foot3", c("Z", "Y")],
            dynamic_t[dynamic_t$marker == "Foot4", c("Z", "Y")]
        ))
    }

    save(x,
        y,
        events,
        ts,
        nL,
        nM,
        nD,
        nS,
        nT,
        dT,
        links,
        sampfreq,
        file = sprintf("%s/processed_data.rdata", outputdir)
    )

    return(list(
        x = x,
        y = y,
        events = events,
        ts = ts,
        nL = nL,
        nM = nM,
        nD = nD,
        nS = nS,
        nT = nT,
        dT = dT,
        sampfreq = sampfreq,
        links = links
    ))
}