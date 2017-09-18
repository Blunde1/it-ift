// 2016 Berent Lunde

/**
* Makes the following available for the TMB user:
* The "cType<Type>" complex AD data type. 
* The arithmetic follows standard arithmetic for complex variables.
* It is defined to work together with the standard TMB data type "Type".
* Constructor: cType<Type> z; defaults to 0+0*i.
*              cType<Type> z((Type)Re,(Type)Im) = Re + i*Im.
* Compound assignements (+=, -=, *=, /=).
* Unary operators (+, -).
* Arithmetic (+, -, *, /).
* Relational and comparison operators (==, !=).
* Standard functions for complex variables (abs, arg, conj, real, imag).
* Exponential functions (exp, log).
* Power functions (pow, sqrt).
* Trigonometric functions (sin, cos, tan, asin, acos, atan, sinh, cosh, tanh).
* 
*/

#ifndef __COMPLEX_HPP_INCLUDED__
#define __COMPLEX_HPP_INCLUDED__

template<class Type>
struct cType{
    
    // MEMBER FUNCTIONS
    Type r, i;
    
    // Constructor
    cType(void) {r=0;i=0;}
    cType(Type r_, Type i_=0) : r(r_), i(i_) {}
    cType(const cType& c) : r(c.r), i(c.i) {}
    
    // Compound assignements
    cType& operator =(const Type& t){
        r = t;
        i = 0;
        return *this;
    }
    cType& operator =(const cType& c){
        r = c.r;
        i = c.i;
        return *this;
    }
    cType& operator +=(const Type& t){
        r = r + t;
        return *this;
    }
    cType& operator +=(const cType& c){
        r = r + c.r;
        i = i + c.i;
        return *this;
    }
    cType& operator -=(const Type& t){
        r = r - t;
        return *this;
    }
    cType& operator -=(const cType& c){
        r = r - c.r;
        i = i - c.i;
        return *this;
    }
    cType& operator *=(const Type& t){
        r = r*t;
        i = i*t;
        return *this;
    }
    cType& operator *=(const cType& c){
        Type tmp_r, tmp_i;
        tmp_r = r*c.r - i*c.i;
        tmp_i = r*c.i + i*c.r;
        r = tmp_r, i=tmp_i;
        return *this;
    }
    cType& operator /=(const Type& t){
        r = r / t;
        i = i / t;
        return *this;
    }
    cType& operator /=(const cType& c){
        Type div = c.r*c.r + c.i*c.i, tmp_r, tmp_i;
        tmp_r = (r*c.r + i*c.i)/div;
        tmp_i = (i*c.r - r*c.i)/div;
        r = tmp_r, i = tmp_i;
        return *this;
    }
};

// NON-MEMBER FUNCTIONS
// Unary operators
template<class Type>
cType<Type> operator+ (const cType<Type>& z){
    cType<Type> c = z;
    return c;
}
template<class Type>
cType<Type> operator- (const cType<Type>& z){
    cType<Type> c = z, zero;
    return zero -= c;
}

// Arithmetic
template<class Type>
cType<Type> operator +(const cType<Type>& c, const Type& t){
    cType<Type> res = c;
    return res+=t;
}
template<class Type>
cType<Type> operator +(const Type& t, const cType<Type>& c){
    cType<Type> res = c;
    return res+=t;
}
template<class Type>
cType<Type> operator +(const cType<Type>& c_1, const cType<Type>& c_2){
    cType<Type> res = c_1;
    return res += c_2;
}
template<class Type>
cType<Type> operator -(const cType<Type>& c, const Type& t){
    cType<Type> res = c;
    return res-=t;
}
template<class Type>
cType<Type> operator -(const Type& t, const cType<Type>& c){
    cType<Type> res(t,0);
    return res -= c;
}
template<class Type>
cType<Type> operator -(const cType<Type>& c_1, const cType<Type>& c_2){
    cType<Type> res = c_1;
    return res -= c_2;
}
template<class Type>
cType<Type> operator *(const cType<Type>& c, const Type& t){
    cType<Type> res = c;
    return res *= t;
}
template<class Type>
cType<Type> operator *(const Type& t, const cType<Type>& c){
    cType<Type> res(t,0);
    return res *= c;
}
template<class Type>
cType<Type> operator *(const cType<Type>& c_1, const cType<Type>& c_2){
    cType<Type> c1 = c_1, c2=c_2;
    return c1 *= c2;
}
template<class Type>
cType<Type> operator /(const cType<Type>& c, const Type& t){
    cType<Type> res = c;
    return res /= t;
}
template<class Type>
cType<Type> operator /(const Type& t, const cType<Type>& c){
    cType<Type> res(t,0);
    return res /= c;
}
template<class Type>
cType<Type> operator /(const cType<Type>& c_1, const cType<Type>& c_2){
    cType<Type> res = c_1;
    return res /= c_2;
}

// Relational and comparison operators
template<class Type>
bool operator ==(const cType<Type>& lhs, const cType<Type>& rhs){
    cType<Type> c_1=lhs, c_2 = rhs;
    if(c_1.r == c_2.r && c_1.i == c_2.i) {return true;}
    else{return false;}
}
template<class Type>
bool operator ==(const cType<Type>& lhs, const Type& rhs){
    cType<Type> c_1=lhs, c_2(rhs,0);
    if(c_1.r == c_2.r && c_1.i == c_2.i) {return true;}
    else{return false;}
}
template<class Type>
bool operator ==(const Type& lhs, const cType<Type>& rhs){
    cType<Type> c_1(lhs,0), c_2=rhs;
    if(c_1.r == c_2.r && c_1.i == c_2.i) {return true;}
    else{return false;}
}
template<class Type>
bool operator !=(const cType<Type>& lhs, const cType<Type>& rhs){
    cType<Type> c_1=lhs, c_2 = rhs;
    return !(c_1==c_2);
}
template<class Type>
bool operator !=(const cType<Type>& lhs, const Type& rhs){
    cType<Type> c_1=lhs, c_2(rhs,0);
    return !(c_1==c_2);
}
template<class Type>
bool operator !=(const Type& lhs, const cType<Type>& rhs){
    cType<Type> c_1(lhs,0), c_2=rhs;
    return !(c_1==c_2);
}

// Standard functions
template<class Type>
Type abs(const cType<Type>& z){
    cType<Type> c = z;
    Type abs = sqrt(c.r*c.r + c.i*c.i);
    return abs;
}
template<class Type>
Type arg(const cType<Type>& z){ // Returns Arg(z)
    cType<Type> c = z;
    Type arg;
    if(c.i==0 && c.r>0){arg = 0; }
    else if(c.i==0 && c.r<0){arg = (Type)M_PI; }
    else{ arg = atan2(c.i,c.r);}
    return arg;
}
template<class Type>
cType<Type> conj(const cType<Type>& z) {
    cType<Type> c = z;
    c.i = - c.i;
    return c;
}
template<class Type>
Type real(const cType<Type>& z){
    cType<Type> c = z;
    return c.r;
}
template<class Type>
Type imag(const cType<Type>& z){
    cType<Type> c = z;
    return c.i;
}

// Exponential functions
template<class Type>
cType<Type> exp(const cType<Type>& z) {
    cType<Type> c = z;
    Type temp = exp(c.r)*cos(c.i);
    c.i = exp(c.r)*sin(c.i);
    c.r = temp;
    return c;
}
template<class Type>
cType<Type> log(const cType<Type>& z){
    cType<Type> c = z;
    Type r = abs(z), theta = arg(z);
    c.r = log(r);
    c.i = theta;
    return c;
}

// Power functions
template<class Type>
cType<Type> pow(const cType<Type>& z1, const cType<Type>& z2){
    cType<Type> c1=z1, c2=z2, c3;
    Type a=c1.r,b=c1.i, c=c2.r,d=c2.i;
    c3.r = pow((a*a+b*b),(c/2))*exp(-d*arg(c1))*
        (cos(c*arg(c1)+0.5*d*log(a*a+b*b)));
    c3.i = pow((a*a+b*b),(c/2))*exp(-d*arg(c1))*
        (sin(c*arg(c1)+0.5*d*log(a*a+b*b)));
    return c3;
}
template<class Type>
cType<Type> pow(const cType<Type>& z, const Type& t){
    cType<Type> c1=z, c2(t,0);
    c1 = pow(c1,c2);
    return c1;
}
template<class Type>
cType<Type> pow(const Type& t, const cType<Type>& z){
    cType<Type> c1=z, c2(t,0);
    c1 = pow(c2,c1);
    return c1;
}
template<class Type>
cType<Type> sqrt(const cType<Type>& z){
    cType<Type> c1=z, c2(0.5,0);
    c1 = pow(c1,c2);
    return c1;
}

// Trigonometric functions
template<class Type>
cType<Type> sin(const cType<Type>& z){
    cType<Type> c = z, i(0,1);
    c = (exp(i*c) - exp((-(Type)1)*i*c)) / ((Type)2 * i);
    return c;
}
template<class Type>
cType<Type> cos(const cType<Type>& z){
    cType<Type> c = z, i(0,1);
    c = (exp(i*c) + exp(c/i)) / ((Type)2);
    return c;
}
template<class Type>
cType<Type> tan(const cType<Type>& z){
    cType<Type> c = z;
    c = sin(c) / cos(c);
    return c;
}
template<class Type>
cType<Type> asin(const cType<Type>& z){
    cType<Type> c = z, i(0,1);
    c = (-(Type)1)*i*log( i*c + sqrt((Type)1 - (c*c)) );
    return c;
}
template<class Type>
cType<Type> acos(const cType<Type>& z){
    cType<Type> c = z;
    c = asin((-(Type)1 ) * c) + (Type)(M_PI/2);
    return c;
}
template<class Type>
cType<Type> atan(const cType<Type>& z){
    cType<Type> c = z, i(0,1);
    c = (-(Type)(0.5))*i*(log((Type)1 - (i*c)) - log((Type)1 + (i*c) ) );
    return c;
}
template<class Type>
cType<Type> sinh(const cType<Type>& z){
    cType<Type> c = z, i(0,1);
    c = sin((-(Type)1)*i*c) * i;
    return c;
}
template<class Type>
cType<Type> cosh(const cType<Type>& z){
    cType<Type> c = z, i(0,1);
    c = cos(c/i);
    return c;
}
template<class Type>
cType<Type> tanh(const cType<Type>& z){
    cType<Type> c = z;
    c = sinh(c) / cosh(c);
    return c;
}

#endif // __COMPLEX_HPP_INCLUDED__