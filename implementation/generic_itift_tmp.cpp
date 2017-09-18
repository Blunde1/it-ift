/**
 * Standalone (except for cType) and relatively minimal working example of it-ift methodology
 * Author: Berent Å. S. Lunde
 * Date: 17.09.2017
 */

#include<TMB.hpp>

/* include cType */
#include "../../cpp_math_library/includes/complex.hpp"

/* differentials diffusion */
template<class Type>
vector<Type> differentials_diff(Type x0, Type dt, vector<Type> p, int process_type){
    
    /**
     * Returns derivatives of order 0-2 w.r.t. x0 for the drift and diffusion 
     * coefficients of a time homogeneous diffusion process.
     * x0: value of diffusion at time t.
     * dt: timestep between observations.
     * p: parameter vector.
     * process_type: integer to specify the diffusion process. 
     */
    
    /* Parameter names */
    Type kappa, alpha, sigma, delta;
    
    /* Differentials w.r.t. x0 */
    Type drift, drift_1, drift_2, diffusion, diffusion_1, diffusion_2;
    
    /* Return vector - holds value of differentials */
    vector<Type> deriv(6);
    
    switch(process_type)
    {
    case 1: 
        /* Geometric Brownian Motion */
        kappa = par[0], sigma = par[1];
        drift = kappa * x0, diffusion = sigma*x0;
        drift_1 = kappa, diffusion_1 = sigma;
        drift_2 = 0, diffusion_2 = 0;
        break;
    case 2:
        /* log GBM */
        kappa = par[0], sigma=par[1];
        drift = kappa - 0.5*sigma*sigma;
        drift_1 = 0, drift_2 = 0;
        diffusion = sigma;
        diffusion_1 = 0, diffusion_2 = 0;
    }
    
    deriv[0]=drift,    deriv[1]=drift_1,    deriv[2]=drift_2;
    deriv[3]=diffusion, deriv[4]=diffusion_1, deriv[5]=diffusion_2;
    return deriv;
    
}

/* differentials jd */

/* log_cf function */


/* 1. log_cf struct */
/* 1. 1. IFT - GH */
/* 1. 2. IFT - GL */

/* 2. log_mgf = cgf struct */
/* 2. inner problem */
/* 2. 1. newton optimiser */
/* 2. 2.  */


/* objective function */
Type objective_function<Type>::operator()(){
    
    DATA_VECTOR(Xt)
    DATA_SCALAR(dt);
    
    PARAMETER_VECTOR(par_diff);
    PARAMETER_VECTOR(par_jump);
    
    return 0;
}