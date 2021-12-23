// -------------- Model 1 - Time independent BIK --------------------------
// mod1.stan defines the time independent BIK model in stan syntax
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

    real<lower=0> ts[nT];               // Observed times

    int<lower=1, upper = nL> lidx[N];   // link idx
    int<lower=1, upper=nM> midx[N];     // marker idx
    int<lower=1, upper=nT> tidx[N];     // time idx

    vector[nD] x[nL, nM];               // Known positions of markers in LCS
    vector[nD] y[N];                    // Observed marker positions
    vector[nD] ls[nL];                  // length of each link


}   

parameters {
    vector[nD]r [nT];                                // Hip joint position over time.
    vector<lower=-pi(), upper=pi()>[nL] theta[nT];   // Segment angles over time 
    vector[nT] lsigma;                               // Measurement noise on log scale
}
transformed parameters {
    vector<lower=0>[nT] sigma;                  // Measurement noise on standard scale
    sigma = exp(lsigma);    
}

model {
    vector[nD] alpha[N];                    // Expected marker position
    vector[nD] J[nL, nT];                   // Joint centers for each link 
    matrix[nD, nD] Gam[nL, nT];             // Rotation Matrix

// -------------------- Priors ---------------------------------------
    // Measurement noise
    lsigma ~ normal(rep_vector(0,nT), 10);              
    
    // prior for r
    r[1:nT, 1] ~ normal(rep_vector(0, nT), 1e6);
    r[1:nT, 2] ~ normal(rep_vector(0, nT), 1e6);

    // Thetas uninformative
    for(l in 1:nL){
        theta[1:nT, l] ~ uniform(rep_vector(-pi(),nT),  rep_vector(pi(), nT));

    }


// -------------------- Kinematic Model ---------------------------------------
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