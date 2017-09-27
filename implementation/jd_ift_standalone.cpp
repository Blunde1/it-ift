#include<TMB.hpp>
#include "includes/itift.hpp"

template<class Type>
vector<Type> differentials(Type x0, Type dt, vector<Type> par, int process, int scheme){
    Type kappa, alpha, sigma, delta;
    Type lambda=0, mu=0, nu=0, k=0;
    Type drift, drift_1, drift_2, diffusion, diffusion_1, diffusion_2;
    
    switch(process)
    {
    case 1:
        /* log GBM */
        kappa = par[0], sigma=par[1];
        drift = kappa - 0.5*sigma*sigma;
        drift_1 = 0, drift_2 = 0;
        diffusion = sigma;
        diffusion_1 = 0, diffusion_2 = 0;
        break;
    }
    vector<Type> deriv(9);
    deriv[0]=drift,    deriv[1]=drift_1,    deriv[2]=drift_2;
    deriv[3]=diffusion, deriv[4]=diffusion_1, deriv[5]=diffusion_2;
    deriv[6]=lambda, deriv[7]=mu, deriv[8]=nu;
    return deriv;
}

// Characteristic function
template<class Type>
cType<Type> logcf(cType<Type> s, Type x0, Type dt, vector<Type> par, int process, int scheme, int jump){
    //using namespace NQL;
    // Differentials
    vector<Type> deriv = differentials((Type)x0, (Type)dt, (vector<Type>)par, (int)process, (int)scheme);
    Type drift, drift_1, drift_2, diffusion, diffusion_1, diffusion_2;
    Type lambda, mu, nu;
    drift = deriv[0], drift_1=deriv[1], drift_2=deriv[2];
    diffusion=deriv[3], diffusion_1=deriv[4], diffusion_2=deriv[5];
    lambda=deriv[6], mu=deriv[7], nu=deriv[8];

    // Characteristic functions
    cType<Type> lcf, lcf_RS;
    Type c_1, c_2, c_3, c_4;
    
    cType<Type> i(0,1);
    
    switch (scheme) {
    case 1: 
        /* Euler-Maruyama */
        lcf = (s*i * (x0 + drift*dt) - pow(s * diffusion, (Type)2)*dt / (Type)2);
        break;
    case 2:
        /* Milstein */
        c_1=diffusion, c_2=0.5*diffusion*diffusion_1, c_3=dt*(drift-c_2);
        lcf_RS = (-c_1*c_1*s*s*dt) / ((Type)2 - ((Type)4)*i*s*c_2*dt  );// -
        //((Type)0.5)*log((Type)1 - ((Type)2)*dt*c_2*s*i);
        lcf = lcf_RS + (i*s*(x0+c_3));
        break;
    case 3:
        /* Scheme 3 */
        c_1 = diffusion + (drift*diffusion_1 + 0.5*pow(diffusion,2)*diffusion_2)*dt;
        c_2 = 0.5*diffusion*diffusion_1;
        c_3 = diffusion*drift_1 - drift*diffusion_1 - 0.5*pow(diffusion, 2)*diffusion_2;
        c_4 = (drift - 0.5*diffusion*diffusion_1)*dt + 
            (drift*drift_1 + 0.5*pow(diffusion, 2)*drift_2)*0.5*pow(dt, 2);
        lcf_RS = -dt*s*s*(((Type)6)*c_1*c_1 + ((Type)6)*c_1*c_3*dt + ((Type)2)*c_3*c_3*dt*dt - i*s*c_3*c_3*dt*dt*dt*c_2) / ((Type)12 - 
            ((Type)24)*i*s*c_2*dt) - ((Type)0.5)*log((Type)1 - ((Type)2)*dt*c_2*s*i);
        lcf = lcf_RS + (i*s*(x0+c_4));
        break;
    }
    
    switch(jump){
    case 0:
        lcf = lcf; // Nothing happens
        break;
    case 1: // Normally distirbuted jumps
        lcf = lcf + lambda*dt*(exp(s*i*mu - ((Type)0.5)*nu*nu*s*s)- (Type)1);
        break;
    case 2: // Gamma distributed jumps
        // mu>0: shape, nu>0: scale
        lcf = lcf + lambda*dt*(pow(((Type)1 - nu*i*s), -mu) - (Type)1);
        break;
    }
    
    return lcf;
}

template<class Type>
struct logcf_f{
    Type x0, dt;
    vector<Type> par;
    int process, scheme, jump;
    logcf_f(Type x0_, Type dt_, vector<Type> par_, int process_, int scheme_, int jump_) :
        x0(x0_), dt(dt_), par(par_), process(process_), scheme(scheme_), jump(jump_) {}
    void set_x0(Type x0_){ this->x0=x0_;}
    template <typename T>
    cType<T> operator()(cType<T> s) { //(const cType<Type>& s) {
        T x0_ = T(x0);
        T dt_ = T(dt);
        vector<T> par_ = par.template cast<T>();
        cType<T> logcf_ = logcf(s, x0_, dt_, par_, (int)process, (int)scheme, (int)jump);
        //cType<T> i(0,1);
        //cType<T> logcf_ = i*s*T(1) - T(0.5)*s*s*T(4); // N(1,2) lcf for checking Gauss-Hermite ift
        return logcf_;
    }
};

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

template<class Type, class Functor>
Type ift_gulag(Type x, Functor f, matrix<Type> quadRules){
    cType<Type> i(0,1);
    Type sum = 0; 
    int n = quadRules.rows();
    vector<Type> ai = quadRules.col(0);
    vector<Type> w = quadRules.col(1);
    for(int j=0; j<n; j++){
        sum += w[j]*real(exp(f(cType<Type>(ai[j])) - i*ai[j]*x+ai[j]));
    }
    sum = sum /((Type)M_PI);
    return sum;
}

template<class Type>
Type objective_function<Type>::operator()()
{
    DATA_VECTOR(x);
    DATA_SCALAR(dt);
    DATA_INTEGER(process);
    DATA_INTEGER(scheme);
    DATA_INTEGER(jump);
    DATA_INTEGER(ghiter);
    DATA_INTEGER(quadrature); // 1:hermite, 2:laguerre
    PARAMETER_VECTOR(par);
    
    Type nll = 0;
    
    // lcf object
    logcf_f<Type> lcf(x[0], dt, par, process, scheme, jump);

    // Qaadrature rules    
    matrix<Type> rules(ghiter,2);
    switch(quadrature){
    case 1: // Hermite
        rules = gauss_hermite_quad(rules);
        for(int k=1; k<x.size(); k++){
            lcf.set_x0(x[k-1]);
            nll -= log(ift_gauher_n(x[k],lcf,rules));
            }
        break;
    case 2: // Laguerre
        rules = gauss_laguerre_quad(rules);
        for(int k=1; k<x.size(); k++){
            lcf.set_x0(x[k-1]);
            nll -= log(ift_gulag(x[k], lcf, rules));
        }
        break;
    }
    
    return nll;
    
}