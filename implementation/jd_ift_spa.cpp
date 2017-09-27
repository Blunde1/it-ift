/**
 * Generic spa methodology - works for diffusions
 * Author: Berent Å. S. Lunde
 * Date: 26.09.2017
 */

#include<TMB.hpp>
#include<fenv.h>
#include "includes/itift.hpp"

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

    // Parameters
    PARAMETER_VECTOR(par); // jump diffusion parameters
    
    int nobs=X.size(), j;
    Type nll = 0;
    vector<Type> s(1);
    //s.setOnes();
    s(0) = Type(10);
    
    // Create inner problem
    spa_iprob<Type> iprob(X(1), X(0), dt, par, process, scheme, jump);
    
    for(j=1; j<nobs; j++){
        iprob.set_x0(X(j-1));
        iprob.set_x( X(j) );
        s(0) = newton_local_extrema(iprob, s, niter);
        nll -= lspa(iprob, s);
    }
    
    return nll;
    
}