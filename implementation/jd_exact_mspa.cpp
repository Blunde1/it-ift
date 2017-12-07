/**
 * Generic mixed exact spa methodology
 * Author: Berent Ã. S. Lunde
 * Date: 07.12.2017
 */

#include<TMB.hpp>
//#include<fenv.h>
#include "includes/itift.hpp"

/* objective function */
template<class Type>
Type objective_function<Type>::operator()(){
    
    // DATA
    DATA_VECTOR(X);
    DATA_SCALAR(dt);
    DATA_INTEGER(process);
    DATA_INTEGER(scheme);
    DATA_INTEGER(jump);
    DATA_INTEGER(niter_no_jump);
    DATA_INTEGER(niter_jump);
    DATA_INTEGER(ghiter_no_jump);
    DATA_INTEGER(ghiter_jump);
    DATA_INTEGER(line_search);
    
    // PARAMETERS
    PARAMETER_VECTOR(par); // jump diffusion parameters
    
    // DECLARATIONS
    int nobs=X.size(), j;
    int par_size = par.size();
    Type nll = Type(0), lfx_no_jump=Type(0), lfx_jump=Type(0), mspa=Type(0), alpha=Type(1);
    vector<Type> s_no_jump(1), s_jump(1), deriv(6);
    s_no_jump.setZero();
    s_jump.setZero();
    
    // Collect lambda
    Type lambda = Type(0);
    switch(jump){
    case(11):
        lambda = exp(par(par_size-3));
        break;
    }
    
    // QUADRATURE RULES
    matrix<Type> rules_no_jump(ghiter_no_jump,2);
    matrix<Type> rules_jump(ghiter_jump,2);
    rules_no_jump = gauss_hermite_quad(rules_no_jump);
    rules_jump = gauss_hermite_quad(rules_jump);
    
    // CGFs
    cgf_s<Type> cgf_no_jump(X(0), dt, par, process, scheme, 0);
    cgf_s<Type> cgf_jump(X(0), dt, par, process, scheme, jump);
    
    // INNER PROBLEMs
    spa_iprob<Type> iprob_no_jump(cgf_no_jump, X(1));
    spa_iprob<Type> iprob_jump(cgf_jump, X(1));
    
    // BUILD NLL
    for(j=1; j<nobs; j++){

        // Update cgf
        cgf_no_jump.set_x0(X(j-1));
        cgf_jump.set_x0(X(j-1));

        // Update inner problem
        iprob_no_jump.set_x0(X(j-1));
        iprob_jump.set_x0(X(j-1));
        iprob_no_jump.set_x( X(j) );
        iprob_jump.set_x( X(j) );
        
        // Solve inner problem
        if(line_search == 0){
            s_no_jump(0) = newton_local_extrema(iprob_no_jump, s_no_jump, niter_no_jump, alpha); // start Newton with 0
            s_jump(0) = newton_local_extrema(iprob_jump, s_jump, niter_jump, alpha);
        }else{
            s_no_jump(0) = line_search_newton(iprob_no_jump, s_no_jump, alpha, niter_no_jump);
            s_no_jump(0) = newton_local_extrema(iprob_no_jump, s_no_jump, 2, alpha); // Do two more Newton.. should not be necessary...
        }
        
        // Create standardized lcf
        lcf_standardized<Type> lphi_01_no_jump(cgf_no_jump, X(j), s_no_jump);
        lcf_standardized<Type> lphi_01_jump(cgf_jump, X(j), s_jump);
        
        // Calculate log exact spa
        lfx_no_jump = log_exact_spa(X(j), cgf_no_jump, s_no_jump, lphi_01_no_jump, rules_no_jump );
        lfx_jump = log_exact_spa(X(j), cgf_jump, s_jump, lphi_01_jump, rules_jump );
        
        // Update likelihood
        mspa = exp(lfx_no_jump - lambda*dt) + exp(lfx_jump)*(Type(1) - exp(-lambda*dt));
        nll -= log(mspa);
    }
    
    return nll;
}