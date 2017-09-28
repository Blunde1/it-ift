/**
 * Generic exact spa methodology
 * Author: Berent Å. S. Lunde
 * Date: 28.09.2017
 */

#include<TMB.hpp>
#include<fenv.h>
#include "includes/itift.hpp"

template<class Type>
struct lcf_standardized{
    spa_iprob<Type> iprob;
    Type x;
    vector<Type> sp; // saddlepoint
    
    lcf_standardized(spa_iprob<Type> iprob_, Type x_, vector<Type> sp_) : 
        iprob(iprob_), x(x_), sp(sp_) {}
    
    template<typename T>
    cType<T> operator()(cType<T> s){
        T x_ = T(x);
        vector<T> sp_ = sp.template cast<T>();
        matrix<T> H = autodiff::hessian(iprob, sp_);
        cType<T> res, i(0,1);
        res = -iprob(sp_)+sp_(0)*x_  - i*s*x_ / sqrt(abs(H(0,0))) + 
            iprob(
                sp_(0) + s*i / sqrt(abs(H(0,0)))
            )
            + (sp_(0) + s*i / sqrt(abs(H(0,0))))*x_;
        return res;
    }
};

/* objective function */
template<class Type>
Type objective_function<Type>::operator()(){
    feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); // Extra line needed
    
    // Data
    DATA_VECTOR(X)
    DATA_SCALAR(dt);
    DATA_INTEGER(process);
    DATA_INTEGER(scheme);
    DATA_INTEGER(jump);
    DATA_INTEGER(niter);
    DATA_INTEGER(qhiter);
    
    // Parameters
    PARAMETER_VECTOR(par); // jump diffusion parameters
    
    int nobs=X.size(), j;
    Type nll = 0, fx, fx_standardized;
    vector<Type> s(1), deriv(6);
    
    // Qaadrature rules    
    matrix<Type> rules(qhiter,2);
    rules = gauss_hermite_quad(rules);
    
    // Create inner problem
    spa_iprob<Type> iprob(X(1), X(0), dt, par, process, scheme, jump);

    // Build nll from spa
    for(j=1; j<nobs; j++){
        // Update inner problem
        iprob.set_x0(X(j-1));
        iprob.set_x( X(j) );
        // Solve inner problem
        deriv = differentials_diff(X(j-1), dt, par, process);
        s(0) = (X(j) - X(j-1) - deriv(0)*dt) / (deriv(3)*deriv(3)*dt); // start at spa for Normal approximation
        s(0) = newton_local_extrema(iprob, s, niter);
        // Create standardized lcf with a get_value method
        lcf_standardized<Type> lphi_01(iprob, X(j), s);
        // Invert to obtain standardized density
        fx_standardized = ift_gauher(Type(0), lphi_01, rules);
        //REPORT(fx_standardized);
        // Calculate exact spa
        fx = iprob(s) + log(fx_standardized) - 
            Type(0.5)*atomic::logdet(autodiff::hessian(iprob,s));
        nll -= log(fx);
    }
    
    return nll;
    
}