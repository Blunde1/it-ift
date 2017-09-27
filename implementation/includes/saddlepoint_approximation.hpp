
#ifndef __SADDLEPOINT_APPROXIMATION_HPP_INCLUDED__
#define __SADDLEPOINT_APPROXIMATION_HPP_INCLUDED__

#include "generating_functions.hpp"
//#include "optimization.hpp"

/* inner problem for saddlepoint approximation - struct */
template<class Type>
struct spa_iprob{
    Type x, x0, dt;
    vector<Type> param;
    int process, scheme, jump;
    spa_iprob(Type x_, Type x0_, Type dt_, vector<Type> param_, 
              int process_, int scheme_, int jump_) : 
        x(x_), x0(x0_), dt(dt_), par(par_), process(process_), 
        scheme(scheme_), jump(jump_) {}
    template<typename T>
    T operator()(vector<T> s){
        T cgf = cgf_fun(s[0], T(x0), T(dt), 
                        par.template cast<T>(),
                        (int)process, (int)scheme, (int)jump);
        return cgf - s[0]*x;
    }
};

/* spa function */
/**
 * returns the log saddlepoint approximation
 * input: saddlepoint, inner_problem, gird point
 */
template<class Type, class Functor>
Type lspa(Functor inner_problem, vector<Type> s) {
    matrix<Type> H = autodiff::hessian(inner_problem, s); // Has dim = 1
    Type lspa = inner_problem(s) - log(sqrt(2*M_PI*exp(atomic::logdet(H))));
    return lspa;
}


#endif //__SADDLEPOINT_APPROXIMATION_HPP_INCLUDED__