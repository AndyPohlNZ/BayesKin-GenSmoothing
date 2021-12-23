// sim-pendulum.stan
// Solves the simple pendulum ODE via an order 4 RK integrator.


functions {

  // Pendulum ODE as a function of theta and theta_dot
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
    int<lower=1> nT;                // Number of timepoints

    vector[2] state0;                 // Inital state
    real t0;                        // Initial time
    real<lower=0> ts[nT];               // Observed times
    vector[2] parms;                  // Parameters [g, Lm]

}   

model {
   
}

generated quantities {
  vector[2] state[nT];
    state[1] = state0;
    state[2:nT] = ode_rk45(pendulum, state0, t0, ts[2:nT], parms);
}
