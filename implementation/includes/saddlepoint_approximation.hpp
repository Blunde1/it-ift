// Copyright (C) 2016-2017 Berent Lunde

#ifndef __SADDLEPOINT_APPROXIMATION_HPP_INCLUDED__
#define __SADDLEPOINT_APPROXIMATION_HPP_INCLUDED__

#include "generating_functions.hpp"
#include "optimization.hpp"

/* inner problem for saddlepoint approximation - struct */
template<class Type>
struct spa_iprob{
    Type x, x0, dt;
    vector<Type> par;
    int process, scheme, jump;
    
    // Constructors
    /* Initialise */
    spa_iprob(Type x_, Type x0_, Type dt_, vector<Type> par_, 
              int process_, int scheme_, int jump_) : 
        x(x_), x0(x0_), dt(dt_), par(par_), process(process_), 
        scheme(scheme_), jump(jump_) {}
    /* Inherit from another inner problem */
    spa_iprob(const spa_iprob& iprob) : 
        x(iprob.x), x0(iprob.x0), dt(iprob.dt), par(iprob.par), 
        process(iprob.process), scheme(iprob.scheme), jump(iprob.jump) {}
    /* Inherit from cgf constructor */
    spa_iprob(const cgf_s<Type> &cg, Type x_) : 
        x(x_), x0(cg.x0), dt(cg.dt), par(cg.par), 
        process(cg.process), scheme(cg.scheme), jump(cg.jump) {}
    
    // Set methods
    void set_x0(Type x0) { this -> x0 = x0; }
    void set_x(Type x_) { this -> x = x_; }
    
    template<typename T>
    T operator()(vector<T> s){
        vector<T> p = par.template cast<T>();
        T cgf = cgf_fun(s(0), T(x0), T(dt), p,
                        (int)process, (int)scheme, (int)jump);
        return cgf - s(0)*T(x);
    }
    template<typename T>
    cType<T> operator()(cType<T> s){
        vector<T> p = par.template cast<T>();
        cType<T> res, i(0,1);
        res = log_cf_fun(-i*s, T(x0), T(dt), p, 
                         (int)process, (int)scheme, (int)jump);
        return res - s*T(x);
    }
};

/* spa function */
/**
 * returns the log saddlepoint approximation
 * input: saddlepoint, inner_problem, gird point
 */
template<class Type, class Functor>
Type lspa(Functor inner_problem, vector<Type> s) {
    matrix<Type> H = autodiff::hessian(inner_problem, s); // Hessian of inner problem equals that of the cgf
    Type lspa = inner_problem(s) - Type(0.5)*(log(Type(2)*M_PI)+atomic::logdet(H));//log(sqrt(Type(2)*M_PI*exp(atomic::logdet(H))));
    return lspa;
}


/* standardized log characteristic function */
template<class Type>
struct lcf_standardized{
    cgf_s<Type>& cgf;
    Type x;
    vector<Type> sp; // saddlepoint
    
    lcf_standardized(cgf_s<Type> &cgf_, Type x_, vector<Type> sp_) : 
        cgf(cgf_), x(x_), sp(sp_) {}
    
    template<typename T>
    cType<T> operator()(cType<T> s){
        T x_ = T(x);
        vector<T> sp_ = sp.template cast<T>();
        matrix<T> H = autodiff::hessian(cgf, sp_);
        cType<T> res, i(0,1);
        res = -cgf(sp_) - i*s*x_ / sqrt(abs(H(0,0))) + 
            cgf( sp_(0) + s*i / sqrt(abs(H(0,0))) );
        return res;
    }
};

/* Exact saddlepoint approximation */
template<class Type>
Type log_exact_spa(Type x, cgf_s<Type> cgf, vector<Type> sp, 
               lcf_standardized<Type> lcf_stand, matrix<Type> rules){
    Type lfx, fx_standardized;
    // Invert to obtain standardized density
    //fx_standardized = ift_gauher(Type(0), lcf_stand, rules);
    fx_standardized = ift_gh_scaled0(lcf_stand, rules);
    
    // Calculate exact spa
    lfx = cgf(sp) - sp(0)*x + log(fx_standardized) - 
        Type(0.5)*atomic::logdet(autodiff::hessian(cgf,sp));
    
    return lfx;
}


#endif //__SADDLEPOINT_APPROXIMATION_HPP_INCLUDED__