// Model 3: 'True' ODE model
// Created as part of a Bayesian Contnuous Inverse Kinmeatics paper.
functions {
    vector pendulum(real t,             // time
                    vector state,       // state[theta, dtheta]
                    vector parms       // parameters [g, Lm]
    ){
                    vector[2] dydt;

                    dydt[1] = state[2]; 
                    dydt[2] =  -parms[1]/parms[2] * sin(state[1]);
                    return dydt;
                    }
}


data {
    int<lower=1> N;                     // Number of data poitns: nM*nS*nT
    int<lower=1> nM;                    // Number of markers
    int<lower=1> nD;                    // Number of dimensions
    int<lower=1> nT;                    // Number of timepoints
    real<lower=0> ts[nT];               // Observed times

    int<lower=1, upper=nM> midx[N];     // marker idx
    int<lower=1, upper=nT> tidx[N];     // time idx

    vector[nD] x[nM];                   // Known positions of markers in LCS
    vector[nD] y[N];                    // Observations  
}   

parameters {
    // State parameters
    real<lower = -pi(), upper = pi()> theta0; //Inital angle
    real dtheta0;                             //Initial angular velocity 

    // Global constants
    real<lower = 0> lm;
    real lsigma;             // measurement noise (m)
}

transformed parameters {
    vector[nD] r[nT];
    for(t in 1:nT){
        r[t] = rep_vector(0, nD);
    }
    
    vector[2] state0;                   // Inital state 
    vector[2] state[nT];                // state over time
    vector[2] parms;                      // Parameters for ODE [g, Lm]

    parms[1] = 9.80665;                 // Gravity
    parms[2] = lm;
    
    state0[1] = theta0;
    state0[2] = dtheta0;
    
    // Solve ODE for state over time.
    state[1] = state0;
    state[2:nT] =  ode_rk45(pendulum, state0, ts[1], ts[2:nT], parms); // Get soln to ode system

    // Transform state to standard parameters for IK model
    vector[nT] theta = to_vector(state[1:nT, 1]);
    vector[nT] dtheta = to_vector(state[1:nT, 2]);

    // Measurement noise
    vector<lower = 0>[nT] sigma = rep_vector(exp(lsigma), nT);
}

model {
    vector[nD] alpha[N];                    // Expected marker position
    matrix[nD, nD] Gam;                     // Rotation Matrix

// -------------------- Priors ---------------------------------------
    // Measurement noise
    lsigma ~ normal(0, 10);  

    // Inital State
    theta0 ~ uniform(-pi(), pi());    // Weakly informative (centered at 22deg but values between 90, 0)
    dtheta0 ~ normal(0, 1e6);       // Uninformative informative on negative angles.
    
    // pendulum length
    lm ~ normal(1, 1.0/3);                  // Weakly informative (centered at 1m with data within +/- 1m 99% of
// -------------------- Kinematic Model ---------------------------------------
    // Marker Model
    for(i in 1:N){
        Gam[1,1] = cos(theta[tidx[i]]); Gam[1,2] = -sin(theta[tidx[i]]);
        Gam[2,1] = sin(theta[tidx[i]]); Gam[2,2] = cos(theta[tidx[i]]);
        alpha[i] = r[tidx[i]] + (Gam*x[midx[i]]);
        y[i] ~ normal(alpha[i], sigma[tidx[i]]);
    }
}
