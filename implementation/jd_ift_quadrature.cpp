
/**
 * Generic it-ift methodology using Gauss-Laguerre or Hermite quadrature
 * Author: Berent Å. S. Lunde
 * Date: 26.09.2017
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
    DATA_INTEGER(qiter);
    DATA_INTEGER(quadrature); // 1:hermite, 2:laguerre
    
    // Parameters
    PARAMETER_VECTOR(par); // jump diffusion parameters

    int nobs=X.size(), j;
    Type nll = 0;
    
    // lcf object
    log_cf<Type> lcf(X[0], dt, par, process, scheme, jump);
    
    // Qaadrature rules    
    matrix<Type> rules(qiter,2);
    switch(quadrature){
    case 1: // Hermite
        rules = gauss_hermite_quad(rules);
        for(j=1; j < nobs; j++){
            lcf.set_x0(X[j-1]);
            nll -= log(ift_gauher(X[j],lcf,rules));
        }
        break;
    case 2: // Laguerre
        rules = gauss_laguerre_quad(rules);
        for(j=1; j < nobs; j++){
            lcf.set_x0(X[j-1]);
            nll -= log(ift_gaulag(X[j], lcf, rules));
        }
        break;
    }
    
    return nll;
   
}