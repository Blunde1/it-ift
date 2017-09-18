
#ifndef __FOURIER_INVERSION_HPP_INCLUDED__
#define __FOURIER_INVERSION_HPP_INCLUDED__

/** Gauss laguerre fourier transformation from logarithmich characteristic function to density function 
 */
template<class Type, class Functor>
Type logcf_to_dens(Type x, Type x0, Functor f, matrix<Type> quadRules){
    cType<Type> i(0,1);
    Type sum = 0; 
    int n = quadRules.rows();
    vector<Type> ai = quadRules.col(0);
    vector<Type> w = quadRules.col(1);
    for(int j=0; j<n; j++){
        sum += w[j]*real(exp(f(cType<Type>(ai[j]),x0) - i*ai[j]*x+ai[j]));
    }
    sum = sum /((Type)M_PI);
    return sum;
}

/** Gauss laguerre fourier transformation from logarithmich characteristic function to distribution function 
 */
template<class Type, class Functor>
Type logcf_to_dist(Type x, Type x0, Functor lcf, matrix<Type> quadRules){
    cType<Type> i(0,1);
    cType<Type> Fx(0);
    int n = quadRules.rows();
    vector<Type> ai = quadRules.col(0);
    vector<Type> w = quadRules.col(1);
    for(int j=0; j<n; j++){
        Fx -= w[j]*(exp(lcf(cType<Type>(ai[j]),x0) - i*x*ai[j] + ai[j] ) - 
            exp(lcf(cType<Type>(-ai[j]),x0) + i*x*ai[j] + ai[j] )) / (i*ai[j]);
    }
    Fx = Fx /((Type)2*M_PI);
    Fx += 0.5;
    return real(Fx);
}

/** Gauss hermite fourier transfomr from lcf to density. lcf takes s, ls0, v0, as argument
 */
template<class Type, class Functor>
Type ift_gauher(Type x, Type x0, Functor f, matrix<Type> rules){
    cType<Type> i(0,1);
    Type sum=0;
    int n = (rules.rows()+1) / 2; // symmetric roots and real part of integrand about origin
    vector<Type> ai = rules.col(0);
    vector<Type> w = rules.col(1);
    for(int j=0; j<n; j++){
        sum += w[j]*real(exp(f(cType<Type>(ai[j]), x0)  - i*ai[j]*x + ai[j]*ai[j] ));
    }
    sum = sum / ((Type)(2*M_PI));
    return sum;
}


#endif //__FOURIER_INVERSION_HPP_INCLUDED__