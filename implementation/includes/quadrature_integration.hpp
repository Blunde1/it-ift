#ifndef __QUADRATURE_INTEGRATION_HPP_INCLUDED__
#define __QUADRATURE_INTEGRATION_HPP_INCLUDED__


/**
 * Numerical Recipes in C - Chapter 4.5
 */

// Gauss-Laguerre quadrature weigts and absciasses
template<class Type>
matrix<Type> gauss_laguerre_quad(matrix<Type> mat){
    int n = mat.rows();
    matrix<Type> resMat(n,2);
    
    Type EPS = 3.0e-50;
    int i, its, j;
    Type ai;
    Type p1, p2, p3, pp, ppp, z, z1;
    vector<Type> x(n), w(n);
    int MAXIT = 20; //max number of iterations
    matrix<Type> returnMat(n,2);
    
    for(i=1; i<=n; i++) {
        if (i == 1){ // Initial guess for smallest root
            z = 1 * 3 / (1+2.4*n); // rest was multiplied with alf == 0;
        }
        else if (i==2){ // Initial guess for second root
            z += 15 / (1 + 2.5*n);
        }
        else {
            ai = (Type)(i-2);
            z += ( (1+2.55*ai) / (1.9*ai) ) * (z - x[i-3]) ; 
        }
        for (its=1; its<=MAXIT; its++) {    // Refinement by NM
            p1 = (Type)1.0;
            p2 = (Type)0.0;
            for (j=1; j<=n; j++){
                p3 = p2;
                p2 = p1;
                p1 = (  (2*j-1-z) * p2 - (j-1)*p3    ) / j;
            }
            pp = (n*p1-n*p2) / z;
            //ppp = ((z-1)*pp - n*p1) / z;
            
            z1 = z;
            // z = z1 - 2*p1*pp / (2*pp*pp -pp*ppp); // Halley's step
            
            z = z1 - p1/pp; // Newton step
            if (abs(z-z1)<= EPS) break;
        }
        x[i-1] = z;
        w[i-1] = -exp(lgamma((Type)n)-lgamma((Type)n))/ (pp*n*p2);
    }
    
    resMat.col(0) = x;
    resMat.col(1) = w;
    
    return resMat;
}

template<class Type>
matrix<Type> gauss_hermite_quad(matrix<Type> mat){
    int n = mat.rows();
    matrix<Type> resMat(n,2);
    
    Type EPS = 3.0e-50;
    int i, its, j, m;
    Type p1, p2, p3, pp, z, z1;
    Type PIM4 = 1/pow(M_PI,1.0/4.0);
    vector<Type> x(n), w(n);
    int MAXIT = 20;

    m = (n+1) / 2; // symmetric roots about origin
    for(i=1; i<=m; i++){
        if(i==1){
            z=sqrt((Type)(2*n+1)) - 1.85575*pow((Type)(2*n+1),-1.6667);
        } else if(i==2){
            z -= 1.14*pow((Type)n,0.426) / z;
        } else if(i==3){
            z = 1.86*z - 0.86*x[0];
        } else if(i==4){
            z = 1.91 * z - 0.91*x[1];
        } else{
            z = 2.0*z-x[i-3];
        }
        for(its=1; its<=MAXIT; its++){
            p1=PIM4;
            p2 = 0.0;
            for(j=1; j<=n; j++){
                p3 = p2;
                p2 = p1;
                p1 = z*sqrt(2.0/j)*p2-sqrt(((Type)(j-1))/j)*p3;
            }
            pp=sqrt((Type)2*n)*p2;
            z1=z;
            z=z1 - p1/pp;
            if(abs(z-z1) <=EPS) break;
        }
        x[i-1] = z;
        x[n-i] = -z;
        w[i-1] = 2.0 / (pp*pp);
        w[n-i] = w[i-1];
    }
    resMat.col(0) = x;
    resMat.col(1) = w;
    
    return resMat;
}




#endif //__QUADRATURE_INTEGRATION_HPP_INCLUDED__