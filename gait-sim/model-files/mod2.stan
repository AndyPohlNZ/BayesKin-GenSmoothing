// -------------- Model 2 - Generalised Smoothing BIK --------------------------
// mod2.stan defines the generalised smoothing BIK model in stan syntax
// Written as part of:
//      A Generalised smoothing appraoch reduces uncertainty in Bayesian Solutions
//      to to time continuous planar inverse kinematics problems.
// 
//      By Andy Pohl
//      Human Perforamnce Labortory 
//      Faculty of Kinesiology
//      University of Calgary
//      December 2021
// ----------------------------------------------------------------- ----------


data {
    int<lower=1> N;                     // Number of data poitns: nL*nM*nS*nT
    int<lower=1> nL;                    // Number of links
    int<lower=1> nM;                    // Number of markers
    int<lower=1> nD;                    // Number of dimensions
    int<lower=1> nT;                    // Number of timepoints
    int<lower = 1> nK;                  // Number of basis functions
    int<lower =1, upper=nK> nK2;        // Number of frequencies
    vector<lower=0>[nT]ts;              // Observed times

    int<lower=1, upper = nL> lidx[N];   // link idx
    int<lower=1, upper=nM> midx[N];     // marker idx
    int<lower=1, upper=nT> tidx[N];     // time idx

    matrix[nT,  nK + 3] bbasis;   // basis functions for rx
    matrix[nT, nK +1] fbasis;     // basis functions for theta/ry
    matrix[nT, nK +1] dfbasis;    // basis functions for dot(theta/ry)

    vector[nD] x[nL, nM];               // Known positions of markers in LCS
    vector[nD] y[N];                    // Observations  
    vector[nD] ls[nL];                  // length of each link

}   

parameters {
    vector[nK+3] beta_rx;                       // Basis coefficients for rx
    real mu_ry;                                 // Mean of ry
    vector[nK] beta_ry;                         // Basis coefficients for ry  
    real<lower=-pi(), upper=pi()> mu_theta[nL]; // mean for theta
    vector[nK] beta_theta[nL];                  // Basis coefficients for thetas
    vector[nK+1] beta_lsigma;                   // Basis coefficients for lsigma   
    real <lower=0> epsilon_lsigma;              // Variance of lsigma
    vector[nT] delta;                           // Random effects for non-centered paramtasiation. - https://mc-stan.org/docs/2_28/stan-users-guide/reparameterization.html
}

transformed parameters {
    vector[nD]r [nT];                           // Hip joint center position over time.
    vector[nL] theta[nT];                       // Segment angles over time.
    vector[nT] mu_lsigma;                       // Scale of measurement noise on log scale
    vector<lower=0>[nT] sigma;                  // Measurement noise on non log scale.

    // Model for measurement noise
    mu_lsigma = fbasis * beta_lsigma;
    sigma = exp(mu_lsigma[1:nT] + delta[1:nT]*epsilon_lsigma); 
    
    // Model for r
    r[1:nT, 1] = to_array_1d(bbasis * beta_rx); // rx
    r[1:nT, 2]= to_array_1d((rep_vector(1, nT) * mu_ry) + (fbasis[1:nT, 2:(nK+1)] * beta_ry)); // ry

    // Model for theta
    for(l in 1:nL){
        theta[1:nT, l] = to_array_1d((rep_vector(1, nT)* mu_theta[l]) + (2*atan(fbasis[1:nT, 2:(nK+1)] * beta_theta[l])));
    }
}

model {
    vector[nD] alpha[N];                    // Expected marker position
    vector[nD] J[nL, nT];                   // Joint venters for each link 
    matrix[nD, nD] Gam[nL, nT];             // Rotation Matrix

// -------------------- Priors ---------------------------------------
    
    // Measurement noise
    epsilon_lsigma ~ uniform(0, 1);
    beta_lsigma ~ normal(rep_vector(0, nK+1), 1e6);   
    delta ~ normal(rep_vector(0,nT), 1);
        
    // prior for r
    beta_rx ~ normal(rep_vector(0, nK+3), 1e6);
    mu_ry ~ normal(0, 1e6);
    beta_ry ~ normal(rep_vector(0,nK), 1e6);

    // Priors for theta spline (uninformative...)
    for(l in 1:nL){
        mu_theta[l] ~ uniform(-pi(), pi());
        beta_theta[l] ~ normal(rep_vector(0,nK), 0.2);  // Uninformative once transformed via 2*tan^{-1}
    }

// -------------------- Inverse Kinematic Model ---------------------------------------
    // Rotation matricies
    for(t in 1:nT){
        for(l in 1:nL){
            Gam[l, t] = to_matrix([cos(theta[t,l]), sin(theta[t,l]), 
            -sin(theta[t,l]), cos(theta[t,l])], nD, nD);
        }
    }

    // Joint centers
    for(t in 1:nT){
        J[1, t] = r[t];
        for(l in 2:nL){
            J[l, t] = J[l-1, t] + (Gam[l-1, t, 1:nD, 1:nD] * ls[l-1]);
        }
    }

    // Marker Model
    for(i in 1:N){
        alpha[i] = J[lidx[i], tidx[i]] + Gam[lidx[i], tidx[i]]*x[lidx[i], midx[i]];
        y[i] ~ normal(alpha[i], sigma[tidx[i]]);
    }
}