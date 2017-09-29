/**
 * Generic exact spa methodology
 * Author: Berent Å. S. Lunde
 * Date: 28.09.2017
 */

#include<TMB.hpp>
#include "includes/itift.hpp"

/* objective function */
template<class Type>
Type objective_function<Type>::operator()(){

    // Data
    DATA_VECTOR(X)
    DATA_SCALAR(dt);
    DATA_INTEGER(process);
    DATA_INTEGER(scheme);
    DATA_INTEGER(jump);
    DATA_INTEGER(niter);
    DATA_INTEGER(ghiter);
    
    // Parameters
    PARAMETER_VECTOR(par); // jump diffusion parameters
    
    int nobs=X.size(), j;
    Type nll = 0, lfx, fx_standardized;
    vector<Type> s(1), deriv(6);
    
    // Qaadrature rules    
    matrix<Type> rules(ghiter,2);
    rules = gauss_hermite_quad(rules);
    
    // Create cgf
    cgf_s<Type> cgf(X(0), dt, par, process, scheme, jump);
    
    // Create inner problem
    //spa_iprob<Type> iprob(X(1), X(0), dt, par, process, scheme, jump);
    spa_iprob<Type> iprob(cgf, X(1));
        
    // Build nll from spa
    for(j=1; j<nobs; j++){
        // Update cgf
        cgf.set_x0(X(j-1));
        // Update inner problem
        iprob.set_x0(X(j-1));
        iprob.set_x( X(j) );
        // Solve inner problem
        deriv = differentials_diff(X(j-1), dt, par, process);
        s(0) = (X(j) - X(j-1) - deriv(0)*dt) / (deriv(3)*deriv(3)*dt); // start at spa for Normal approximation
        s(0) = newton_local_extrema(iprob, s, niter);
        // Create standardized lcf
        lcf_standardized<Type> lphi_01(cgf, X(j), s);
        // Calculate log exact spa
        lfx = log_exact_spa(X(j), cgf, s, lphi_01, rules );
        // Update likelihood
        nll -= lfx;
    }
    
    return nll;
    
}