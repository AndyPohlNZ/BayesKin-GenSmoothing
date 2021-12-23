// obs-pendulum.stan
// Generates marker observations from a simple pendulum.



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
    int<lower=1> nT;            // Number of timepoints
    int<lower=1> nM;            // Number of markers
    int<lower=1> nD;            // Number of dimensions

    vector[2] state0;                   // Inital state 
    real t0;                    // Initial time
    real<lower=0> ts[nT];               // Observed times
    vector[2] parms;              // Parameters

    vector[nD] x[nM];           // Known positions of markers in LCS
    vector[nD] r[nT];           // position of pendulum base point over time
    real<lower=0> sigma;        // measurement noise im meters

}       


model {
   
}

generated quantities {
    vector[2] state[nT];                // state over time

    matrix[nD, nD] Gam;                         // Rotation Matrix
    vector[nD] alpha[nM];                       // Expected marker position
    real alpha_all[nM, nD, nT];             // Storage for alpha over time
    vector[nT] theta;                         // Theta over time.
    real y[nM, nD, nT];                     // y over time.

    // Solve ODE via Runge-Kutta 4th order
    state[1] = state0;
    state[2:nT] = ode_rk45(pendulum, state0, t0, ts[2:nT], parms);

    theta = to_vector(state[1:nT, 1]);  
    
    // Rigid Body model:
    for(t in 1:nT){
        // Construct rotation matrix at time t
        Gam[1,1] = cos(theta[t]); Gam[1,2] = -sin(theta[t]);
        Gam[2,1] = sin(theta[t]); Gam[2,2] = cos(theta[t]);

        for(j in 1:nM){
            // Expected marker position
            alpha[j] = r[t] + (Gam * x[j]);
            alpha_all[j, 1:nD, t] = to_array_1d(alpha[j, 1:nD]);
        
            // Marker Observations at time t.
            y[j, 1, t] = alpha[j,1] + normal_rng(0, sigma);
            y[j, 2, t] = alpha[j,2] + normal_rng(0, sigma);
        }
    }
}