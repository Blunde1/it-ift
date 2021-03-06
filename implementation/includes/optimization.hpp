// Copyright (C) 2016-2017 Berent Lunde

#ifndef __OPTIMIZATION_HPP_INCLUDED__
#define __OPTIMIZATION_HPP_INCLUDED__

template<class Type, class Functor>
Type newton_local_extrema(Functor f, vector<Type> s, int niter, Type alpha) {
    //Min/maximize f
    vector<Type> g, s_old=s, s_new(1);
    matrix<Type> H;
    //Type EPS = 1.0e-3;
    for (int i = 0; i<niter; i++) {
        g = autodiff::gradient(f, s_old); 
        H = autodiff::hessian(f, s_old); 
        s_new -= alpha*(atomic::matinv(H) * g);
        //if(abs(s_new(0)-s_old(0)) <= EPS) break;
        //else s_old = s_new;
        s_old = s_new;
    }
    return s_new(0);
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

// Line search Newton method
/**
 * Tested to work for double, should be careful with AD types (but it passed the test).
 */
template<class Type, class Functor>
Type line_search_newton(Functor f, vector<Type> start, Type& alpha, int maxniter){
    
    int i;
    Type EPS = 3.0e-50, f_eps = 1.0e50;
    vector<Type> g, s_old=start, s_new(1);
    matrix<Type> H;
    
    for(i=0; i<maxniter; i++){
        
        g = autodiff::gradient(f,s_old);
        H = autodiff::hessian(f,s_old);
        s_new -= alpha*(atomic::matinv(H)*g);
        
        if (abs(s_new(0)-s_old(0))<= EPS){
            break;
        } else if(f(s_new) <= (f(s_old) + f_eps) ){
            s_old = s_new;
        } else{ // reset with smaller alpha
            alpha = 0.9*alpha;
            s_old = start;
            i = 0;
        }
    }
    return  s_new(0);
}


#endif //__OPTIMIZATION_HPP_INCLUDED__