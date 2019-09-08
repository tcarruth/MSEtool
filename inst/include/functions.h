#include "ns/ns_cDD.h"
#include "ns/ns_SCA.h"
#include "ns/ns_VPA.h"
#include "ns/ns_SRA_scope.h"

//posfun from ADMB
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * (x - eps) * (x - eps));
  return ans;
}

//posfun2 from ADMB
template<class Type>
Type posfun2(Type x, Type eps, Type &penalty) {
  Type x_eps = x/eps;
  Type denom = 1.;
  denom += x_eps;
  denom += pow(x_eps, 2);
  denom += pow(x_eps, 3);
  Type pen2 = pow(denom, -1);
  pen2 *= eps;
  Type ans = CppAD::CondExpGe(x, eps, x, pen2);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * pow(x - eps, 3));
  return ans;
}


// Calculates analytical solution of a lognormal variable
template<class Type>
Type calc_sigma(vector<Type> I_y, vector<Type> Ipred_y) {
  Type sum_square = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.size();y++) {
    if(!R_IsNA(asDouble(I_y(y))) && I_y(y)>0) {
      sum_square += (log(I_y(y)/Ipred_y(y))) * (log(I_y(y)/Ipred_y(y)));
      n_y += 1.;
    }
  }
  Type sigma = pow(sum_square/n_y, 0.5);
  return sigma;
}

// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(vector<Type> I_y, vector<Type> B_y) {
  Type num = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.size();y++) {
    if(!R_IsNA(asDouble(I_y(y))) && I_y(y)>0) {
      num += log(I_y(y)/B_y(y));
      n_y += 1.;
    }
  }
  Type q = exp(num/n_y);
  return q;
}



//////////// Functions for cDD.h, DD.h, SCA.h
template<class Type>
Type BH_SR(Type SSB, Type h, Type R0, Type SSB0) {
  Type Rpred = 4 * h * R0 * SSB;
  Type den = SSB0 * (1-h);
  den += (5*h-1) * SSB;
  Rpred /= den;
  return Rpred;
}

template<class Type>
Type Ricker_SR(Type SSB, Type h, Type R0, Type SSB0) {
  Type SSBPR0 = SSB0/R0;
  Type expon = 1;
  expon -= SSB/SSB0;
  expon *= 1.25;

  Type Rpred = pow(5*h, expon);
  Rpred *= SSB;
  Rpred /= SSBPR0;
  return Rpred;
}

// Sub-annual time steps for SP and iteratively find the F that predicts the observed catch
template<class Type>
Type SP_F(Type U_start, Type C_hist, Type MSY, Type K, Type n, Type nterm, Type dt, int nstep, int nitF,
          vector<Type> &Cpred, vector<Type> &B, int y, Type &penalty) {

  Type F = -log(1 - U_start);
  if(dt < 1) {
    for(int i=0;i<nitF;i++) {
      Type Catch = 0;
      Type B_next = B(y);

      for(int seas=0;seas<nstep;seas++) {
        Type SP = CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B_next / K * log(B_next/K),
                                   nterm/(n-1) * MSY * (B_next/K - pow(B_next/K, n))) - F * B_next;
        SP *= dt;
        Catch += F * B_next * dt;
        B_next += SP;
      }

      if(i==nitF-1) {
        B(y+1) = B_next;
        Cpred(y) = Catch;
      } else {
        F *= C_hist/Catch;
      }
    }
  } else {
    Cpred(y) = C_hist;
    F = CppAD::CondExpLt(1 - C_hist/B(y), Type(0.025), 1 - posfun(1 - C_hist/B(y), Type(0.025), penalty), C_hist/B(y));
    B(y+1) = B(y) + CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B(y) / K * log(B(y)/K),
      nterm/(n-1) * MSY * (B(y)/K - pow(B(y)/K, n))) - F * B(y);
  }

  return F;
}

