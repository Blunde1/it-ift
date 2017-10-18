/**
 * Generic exact spa methodology
 * Author: Berent Å. S. Lunde
 * Date: 28.09.2017
 */

#include<TMB.hpp>
//#include<fenv.h>
#include "includes/itift.hpp"

/* objective function */
template<class Type>
Type objective_function<Type>::operator()(){
    //feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); // Extra line needed
    // Data
    DATA_VECTOR(X);
    DATA_SCALAR(dt);
    DATA_INTEGER(process);
    DATA_INTEGER(scheme);
    DATA_INTEGER(jump);
    DATA_INTEGER(niter);
    DATA_INTEGER(ghiter);
    DATA_INTEGER(line_search);
    
    // Parameters
    PARAMETER_VECTOR(par); // jump diffusion parameters
    
    int nobs=X.size(), j;
    Type nll = 0, lfx=0, alpha=Type(1);
    vector<Type> s(1), deriv(6);
    s.setZero();
    
    // Qaadrature rules    
    matrix<Type> rules(ghiter,2);
    rules = gauss_hermite_quad(rules);
    
    // Create cgf
    cgf_s<Type> cgf(X(0), dt, par, process, scheme, jump);
    
    // Create inner problem
    spa_iprob<Type> iprob(cgf, X(1));
    
    // Build nll from spa
    for(j=1; j<nobs; j++){
        // Update cgf
        cgf.set_x0(X(j-1));
        // Update inner problem
        iprob.set_x0(X(j-1));
        iprob.set_x( X(j) );
        // Solve inner problem
        //deriv = differentials_diff(X(j-1), dt, par, process);
        //s(0) = (X(j) - X(j-1) - deriv(0)*dt) / (deriv(3)*deriv(3)*dt); // start at spa for Normal approximation
        if(line_search == 0){
            s(0) = newton_local_extrema(iprob, s, niter, alpha); // start Newton with 0
        }else{
            s(0) = line_search_newton(iprob, s, alpha, niter);
            s(0) = newton_local_extrema(iprob, s, 2, alpha); // Do two more Newton.. should not be necessary...
        }
        // Create standardized lcf
        lcf_standardized<Type> lphi_01(cgf, X(j), s);
        // Calculate log exact spa
        lfx = log_exact_spa(X(j), cgf, s, lphi_01, rules );
        
        // Update likelihood
        nll -= lfx;
    }

    return nll;
    
}