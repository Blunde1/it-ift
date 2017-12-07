// Copyright (C) 2016-2017 Berent Lunde

#ifndef __differentials_hpp_included__
#define __differentials_hpp_included__

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
    Type lambda, mu, nu, k;
    
    /* Differentials w.r.t. x0 */
    Type drift, drift_1, drift_2, diffusion, diffusion_1, diffusion_2;
    
    /* Return vector - holds value of differentials */
    vector<Type> deriv(6);
    
    switch(process_type)
    {
    case 1: 
        /* Geometric Brownian Motion */
        kappa = p(0), sigma = p(1);
        drift = kappa * x0, diffusion = sigma*x0;
        drift_1 = kappa, diffusion_1 = sigma;
        drift_2 = Type(0), diffusion_2 = Type(0);
        break;
    case 2:
        /* log GBM */
        kappa = p(0), sigma=p(1);
        drift = kappa - Type(0.5)*sigma*sigma;
        drift_1 = Type(0), drift_2 = Type(0);
        diffusion = sigma;
        diffusion_1 = Type(0), diffusion_2 = Type(0);
        break;
    case 3:
        /* log Merton jump-diffusion */
        kappa=p(0), sigma=exp(p(1)), lambda=exp(p(2)), mu=p(3), nu=exp(p(4));
        k = exp(mu+Type(0.5)*nu*nu) - Type(1);
        drift = kappa - lambda*k - Type(0.5)*sigma*sigma;
        diffusion = sigma;
        drift_1 = Type(0), drift_2 = Type(0);
        diffusion_1 = Type(0), diffusion_2 = Type(0);
        break;
    }
    
    deriv(0)=drift,    deriv(1)=drift_1,    deriv(2)=drift_2;
    deriv(3)=diffusion, deriv(4)=diffusion_1, deriv(5)=diffusion_2;
    return deriv;
    
}

#endif // __differentials_hpp_included__