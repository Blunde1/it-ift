#include<TMB.hpp>
#include "includes/itift.hpp"

// normal lcf fun
template<class Type>
cType<Type> lcf_norm(cType<Type>s, Type mu, Type sigma){
    cType<Type> i(0,1), lcf(0,0);
    lcf = i*s*mu - Type(0.5)*s*s*sigma*sigma;
    return lcf;
}

// normal lcf struct
template<class Type>
struct lcf_norm_f{
    Type mu, sigma;
    lcf_norm_f(Type mu_, Type sigma_) : mu(mu_), sigma(sigma_) {}
    template<typename T>
    cType<T> operator()(cType<T> s){
        T mu_ = T(mu);
        T sigma_ = T(sigma);
        cType<T> lcf = lcf_norm(s, mu_, sigma_);
        return lcf;
    }
};

// gh ift
template<class Type, class Functor>
Type ift_gauher_n(Type x, Functor f, matrix<Type> rules){
    cType<Type> i(0,1);
    Type sum=0;
    int n = (rules.rows()+1) / 2; // symmetric roots and real part of integrand about origin
    vector<Type> ai = rules.col(0);
    vector<Type> w = rules.col(1);
    for(int j=0; j<n; j++){
        sum += w[j]*real(exp(f(cType<Type>(ai[j]))  - i*ai[j]*x + ai[j]*ai[j] ));
    }
    sum = sum / ((Type)(M_PI));
    return sum;
}

// gl ift

// obj fun
template<class Type>
Type objective_function<Type>::operator()(){
    DATA_SCALAR(x);
    DATA_INTEGER(gqiter)
    PARAMETER(mu);
    PARAMETER(sigma);
    
    // check that lcf_norm report correct values
    cType<Type> s(1,1), lcf_test;
    lcf_test = lcf_norm(s,mu,sigma);
    REPORT(lcf_test.r);
    REPORT(lcf_test.i);
    // ok
    
    // Create lcf object
    lcf_norm_f<Type> lcf(mu, sigma);
    Type tmp1, tmp2;
    tmp1 = lcf(s).r;
    tmp2 = lcf(s).i;
    REPORT(tmp1); // reports correct values
    REPORT(tmp2);
    
    // Gauss-Hermite quadrature rules
    Type fx;
    matrix<Type> rules(gqiter,2);
    rules = gauss_hermite_quad(rules);
    
    // Invert lcf object
    fx = ift_gauher_n(x, lcf, rules);
    
    // Test passing lcf function to ift gh via pointers
    cType<Type> (*lcf_pointer)(cType<Type>, Type, Type) = lcf_norm;
    
    return fx;
}