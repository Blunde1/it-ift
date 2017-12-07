// Copyright (C) 2016-2017 Berent Lunde

#ifndef __generating_functions_hpp_included__
#define __generating_functions_hpp_included__

#include "differentials.hpp"
#include "complex.hpp"

/* log_cf function */
template<class Type>
cType<Type> log_cf_fun(cType<Type> s, Type x0, Type dt, vector<Type> p, int process_type, int scheme, int jump){
    
    // Differentials
    vector<Type> deriv = differentials_diff(x0, dt, p, process_type);
    Type drift, drift_1, drift_2, diffusion, diffusion_1, diffusion_2;
    drift = deriv(0), drift_1=deriv(1), drift_2=deriv(2);
    diffusion=deriv(3), diffusion_1=deriv(4), diffusion_2=deriv(5);
    
    // Jump parameters
    Type lambda, mu, nu;
    
    //Characteristic function of schemes from - PW-2012 
    cType<Type> lcf(0,0), lcf_RS(0,0), i(0,1);
    Type c_1, c_2, c_3, c_4;
    switch (scheme) {
    case 1: 
        /* Euler-Maruyama */
        lcf = s*i * (x0 + drift*dt) - ((Type)0.5) * s*s*diffusion*diffusion*dt;
        break;
    case 2:
        /* Milstein */
        c_1=diffusion, c_2=0.5*diffusion*diffusion_1, c_3=dt*(drift-c_2);
        lcf_RS = (-c_1*c_1*s*s*dt) / ((Type)2 - ((Type)4)*i*s*c_2*dt  ) -
            ((Type)0.5)*log((Type)1 - ((Type)2)*dt*c_2*s*i);
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
    }
    
    // Adding jumps
    int p_size = p.size();
    switch(jump){
    case 0:
        lcf = lcf; // Nothing happens
        break;
    case 1: // Normally distirbuted jumps
        lambda = exp(p(p_size-3)), mu = p(p_size-2), nu = exp(p(p_size-1));
        lcf = lcf + lambda*dt*(exp(s*i*mu - ((Type)0.5)*nu*nu*s*s)- (Type)1);
        break;
    case 2: // Gamma distributed jumps
        lambda = p[p_size-3], mu = p[p_size-2], nu = p[p_size-1];
        // mu>0: shape, nu>0: scale
        lcf = lcf + lambda*dt*(pow(((Type)1 - nu*i*s), -mu) - (Type)1);
        break;
    case 11: // Zero-truncated compounded Poisson - Normally distributed jumps
        lambda = exp(p(p_size-3)), mu = p(p_size-2), nu = exp(p(p_size-1));
        lcf = lcf + lambda * dt * ( exp( s*i*mu - Type(0.5)*nu*nu*s*s) - Type(1) + 
            log(Type(1) - exp(-lambda*dt*exp(s*i*mu - Type(0.5)*nu*nu*s*s) )) - 
            log(Type(1) - exp(-lambda*dt)) );
    }
    
    return lcf;
}

/* log_cf struct */
template<class Type>
struct log_cf{
    
    Type x0, dt;
    vector<Type> p, p_jump;
    int process, scheme, jump;
    
    log_cf(Type x0_, Type dt_,  vector<Type> p_, int process_, int scheme_, int jump_) :
        x0(x0_), dt(dt_), p(p_), process(process_), scheme(scheme_), jump(jump_) {}
    
    void set_x0(Type x0) { this -> x0 = x0; }
    
    template<typename T>
    cType<T> operator()(cType<T> s) {
        T dt_ = T(dt);
        T x0_ = T(x0);
        vector<T> p_ = p.template cast<T>();
        cType<T> lcf = log_cf_fun(s, x0_, dt_, p_, (int)process, (int)scheme, (int)jump);
        return lcf;
    }
};

/* CGF function */
template<class Type>
Type cgf_fun(Type s, Type x0, Type dt, vector<Type> p, int process_type, int scheme, int jump){
    cType<Type> i(0,1);
    Type cgf = log_cf_fun(-i*s, x0, dt, p, process_type, scheme, jump).r;
    return cgf;
}

/* cgf struct */
template<class Type>
struct cgf_s{
    
    // Member
    Type x0, dt;
    vector<Type> par;
    int process, scheme, jump;
    
    // Constructors
    cgf_s(Type x0_, Type dt_, vector<Type> par_, int process_, int scheme_, int jump_) : 
        x0(x0_), dt(dt_), par(par_), process(process_), scheme(scheme_), jump(jump_) {}
    cgf_s(const cgf_s& cf) : 
        x0(cf.x0), dt(cf.dt), par(cf.par), process(cf.process), scheme(cf.scheme), jump(cf.jump) {}
    
    // Set methods
    void set_x0(Type x0) { this -> x0 = x0; }
    
    // Operators
    template<typename T>
    T operator()(T s){
        T dt_ = T(dt);
        T x0_ = T(x0);
        vector<T> p_ = par.template cast<T>();
        T cgf = cgf_fun(s, x0_, dt_, p_, (int)process, (int)scheme, (int)jump);
        return cgf;
    }
    template<typename T>
    T operator()(vector<T> s){
        T dt_ = T(dt);
        T x0_ = T(x0);
        vector<T> p_ = par.template cast<T>();
        T cgf = cgf_fun(s(0), x0_, dt_, p_, (int)process, (int)scheme, (int)jump);
        return cgf;
    }
    template<typename T>
    cType<T> operator()(cType<T> s){
        T dt_ = T(dt);
        T x0_ = T(x0);
        vector<T> p_ = par.template cast<T>();
        cType<T> res, i(0,1);
        res = log_cf_fun(-i*s, x0_, dt_, p_, 
                         (int)process, (int)scheme, (int)jump);
        return res;
    }
    
};


#endif //__generating_functions_hpp_included__