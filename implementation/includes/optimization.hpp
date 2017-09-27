// Copyright (C) 2016-2017 Berent Lunde

#ifndef __OPTIMIZATION_HPP_INCLUDED__
#define __OPTIMIZATION_HPP_INCLUDED__

template<class Type, class Functor>
vector<Type> newton_local_extrema(Functor f, vector<Type> s, int niter) {
    //Min/maximize f
    for (int i = 0; i<niter; i++) {
        vector<Type> g = autodiff::gradient(f, s); 
        matrix<Type> H = autodiff::hessian(f, s); 
        s -= atomic::matinv(H) * g;
    }
    return s;
}

template<class Type, class Functor>
Type newton_zero(Functor f, vector<Type> &s, int niter) {
    //Find zero of f(s)
    for (int i = 0; i<niter; i++) {
        vector<Type> g = autodiff::gradient(f, s); 
        s[0] -= f(s)/g[0];
    }
    return s[0];
} 


#endif //__OPTIMIZATION_HPP_INCLUDED__