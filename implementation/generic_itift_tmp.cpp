/**
 * Standalone (except for cType) and relatively minimal working example of it-ift methodology
 * Author: Berent Å. S. Lunde
 * Date: 17.09.2017
 */

#include<TMB.hpp>

/* include cType */
#include "includes/itift.hpp"

/* differentials diffusion */
template<class Type>
vector<Type> differentials_diff(Type x0, Type dt, vector<Type> p, vector<Type> p_jump, int process_type){
    
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
        kappa = p[0], sigma = p[1];
        drift = kappa * x0, diffusion = sigma*x0;
        drift_1 = kappa, diffusion_1 = sigma;
        drift_2 = 0, diffusion_2 = 0;
        break;
    case 2:
        /* log GBM */
        kappa = p[0], sigma=p[1];
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
template<class Type>
cType<Type> log_cf_fun(cType<Type> s, Type x0, Type dt, vector<Type> p, vector<Type> p_jump, int process_type, int scheme, int jump){

    // Differentials
    vector<Type> deriv = differentials_diff(x0, dt, p, p_jump, process_type);
    Type drift, drift_1, drift_2, diffusion, diffusion_1, diffusion_2;
    drift = deriv[0], drift_1=deriv[1], drift_2=deriv[2];
    diffusion=deriv[3], diffusion_1=deriv[4], diffusion_2=deriv[5];
    
    // Jump parameters
    Type lambda, mu, nu;
    
    //Characteristic function of schemes from - PW-2012 
    cType<Type> lcf, lcf_RS;
    Type c_1, c_2, c_3, c_4;
    cType<Type> i(0,1);
    switch (scheme) {
    case 1: 
        /* Euler-Maruyama */
        lcf = (s*i * (x0 + drift*dt) - pow(s * diffusion, (Type)2)*dt / (Type)2);
        break;
    case 2:
        /* Milstein */
        c_1=diffusion, c_2=0.5*diffusion*diffusion_1, c_3=dt*(drift-c_2);
        lcf_RS = (-c_1*c_1*s*s*dt) / ((Type)2 - ((Type)4)*i*s*c_2*dt  );// -
        //((Type)0.5)*log((Type)1 - ((Type)2)*dt*c_2*s*i);
        lcf = lcf_RS + (i*s*(x0+c_3));
        break;
    case 3:
        /* Scheme 3 */
        c_1 = diffusion + (drift*diffusion_1 + 0.5*pow(diffusion,2)*diffusion_2)*dt;
        c_2 = 0.5*diffusion*diffusion_1;
        c_3 = diffusion*drift_1 - drift*diffusion_1 - 0.5*pow(diffusion, 2)*diffusion_2;
        c_4 = (drift - 0.5*diffusion*diffusion_1)*dt + 
            (drift*drift_1 + 0.5*pow(diffusion, 2)*drift_2)*0.5*pow(dt, 2);
        lcf_RS = -dt*s*s*(((Type)6)*c_1*c_1 + ((Type)6)*c_1*c_3*dt + ((Type)2)*c_3*c_3*dt*dt - i*s*c_3*c_3*dt*dt*dt*c_2) / ((Type)12 - 
            ((Type)24)*i*s*c_2*dt) - ((Type)0.5)*log((Type)1 - ((Type)2)*dt*c_2*s*i);
        lcf = lcf_RS + (i*s*(x0+c_4));
        break;
    }
    
    // Adding jumps
    switch(jump){
    case 0:
        lcf = lcf; // Nothing happens
        break;
    case 1: // Normally distirbuted jumps
        lambda = p_jump[0], mu = p_jump[1], nu = p_jump[2];
        lcf = lcf + lambda*dt*(exp(s*i*mu - ((Type)0.5)*nu*nu*s*s)- (Type)1);
        break;
    case 2: // Gamma distributed jumps
        lambda = p_jump[0], mu = p_jump[1], nu = p_jump[2];
        // mu>0: shape, nu>0: scale
        lcf = lcf + lambda*dt*(pow(((Type)1 - nu*i*s), -mu) - (Type)1);
        break;
    }
    
    return lcf;
}

/* 1. log_cf struct */
/* 1. 1. IFT - GH */
/* 1. 2. IFT - GL */
template<class Type>
struct log_cf{
    //cType<Type>& s;
    Type x0, dt;
    vector<Type> p, p_jump;
    int process, scheme, jump;
    logcf_f(Type x0_, Type dt_, vector<Type> p_, vector<Type> p_jump_, int process_, int scheme_, int jump_) :
        x0(x0_), dt(dt_), p(p_), p_jump(p_jump_), process(process_), scheme(scheme_), jump(jump_) {}
    template <class T>
    cType<T> operator()(const cType<T>& s, const T& x0) {
        T dt_ = T(dt);
        vector<T> p_ = p.template cast<T>();
        vector<T> p_jump_ = p.template cast<T>();
        cType<T> lcf = log_cf_fun(s, x0, dt_, p_, p_jump_, (int)process, (int)scheme, (int)jump);
        //cType<Type> lcf = log_cf_fun(s, x0_, (Type)dt, (vector<Type>)p, (vector<Type>)p_jd (int)process, (int)scheme, (int)jump);
        return lcf;
    }
};

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