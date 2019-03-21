//posfun from ADMB
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * pow(x - eps, 2));
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
      sum_square += pow(log(I_y(y)/Ipred_y(y)), 2);
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


//////////// Functions for SCA.h
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

template<class Type>
vector<Type> calc_NPR(Type F, vector<Type> vul, vector<Type> M, int max_age) {
  vector<Type> NPR(max_age);
  NPR(0) = 1.;
  for(int a=1;a<max_age;a++) NPR(a) = NPR(a-1) * exp(-vul(a-1) * F - M(a-1));
  NPR(max_age-1) /= 1 - exp(-vul(max_age-1) * F - M(max_age-1)); // Plus-group
  return NPR;
}


template<class Type>
Type sum_EPR(vector<Type> NPR, vector<Type> weight, vector<Type> mat) {
  Type answer = 0.;
  for(int a=0;a<NPR.size();a++) answer += NPR(a) * weight(a) * mat(a);
  return answer;
}

template<class Type>
Type sum_BPR(vector<Type> NPR, vector<Type> weight) {
  Type answer = 0.;
  for(int a=0;a<NPR.size();a++) answer += NPR(a) * weight(a);
  return answer;
}

template<class Type>
Type sum_VBPR(vector<Type> NPR, vector<Type> weight, vector<Type> vul) {
  Type answer = 0.;
  for(int a=0;a<NPR.size();a++) answer += NPR(a) * weight(a) * vul(a);
  return answer;
}

template<class Type>
vector<Type> calc_logistic_vul(vector<Type> vul_par, int max_age, Type &prior) {
  vector<Type> vul(max_age);
  Type maxage = max_age;
  Type a_95 = invlogit(vul_par(0)) * 0.75 * maxage;
  Type a_50 = a_95 - exp(vul_par(1));

  prior -= dnorm(vul_par(1), Type(0), Type(3), true);

  for(int a=0;a<max_age;a++) {
    Type aa = a;
    aa += 1;
    vul(a) = pow(1 + exp(-log(Type(19.0)) * (aa - a_50)/(a_95 - a_50)), -1);
  }
  Type interim_vmax = max(vul);
  for(int a=0;a<max_age;a++) vul(a) /= interim_vmax;

  return vul;
}

template<class Type>
Type dnorm_vul(Type x, Type mu, Type sd) {
  Type res = -0.5;
  res *= pow(x - mu, 2);
  res /= pow(sd, 2);
  return exp(res);
}

template<class Type>
vector<Type> calc_dome_vul(vector<Type> vul_par, int max_age, Type &prior) {
  vector<Type> vul(max_age);

  Type maxage = max_age;
  Type a_full = invlogit(vul_par(0)) * 0.75 * maxage;
  Type a_50 = a_full - exp(vul_par(1));
  Type a_full2 = invlogit(vul_par(2));
  a_full2 *= maxage - a_full;
  a_full2 += a_full;
  Type vul_max = invlogit(vul_par(3));

  //prior -= dnorm(vul_par(0), Type(0), Type(3), true);
  prior -= dnorm(vul_par(1), Type(0), Type(3), true);
  prior -= dnorm(vul_par(2), Type(0), Type(3), true);
  //prior -= dnorm(vul_par(3), Type(0), Type(3), true);

  Type var_asc = pow(a_50 - a_full, 2);
  var_asc /= log(Type(4));

  Type var_des = pow(maxage - a_full2, 2);
  var_des /= -2 * log(vul_max);

  Type sd_asc = pow(var_asc, 0.5);
  Type sd_des = pow(var_des, 0.5);

  for(int a=0;a<max_age;a++) {
    Type aa = a;
    aa += 1;
    Type vul_asc = dnorm_vul(aa, a_full, sd_asc);
    Type vul_des = dnorm_vul(aa, a_full2, sd_des);

    vul(a) = CppAD::CondExpLe(aa, a_full, vul_asc, CppAD::CondExpLe(aa, a_full2, Type(1), vul_des));
  }
  Type interim_vmax = max(vul);
  for(int a=0;a<max_age;a++) vul(a) /= interim_vmax;


  return vul;
}

template<class Type>
Type dlnorm_comp(vector<Type> obs, vector<Type> pred) {
  Type log_lik = 0.;
  for(int a=0;a<obs.size();a++) log_lik += dnorm(log(obs(a)), log(pred(a)), pow(0.01/pred(a), 0.5), true);
  return log_lik;
}


// Iterative solver for F in VPA model -
// if F_ya = Z_ya C_ya / (exp(Z_ya) - 1) N_y+1,a+1, then:
// g(x) = log(Z) + log(C) - log(exp(Z) - 1) - log(N) - x
// where x = log(F), and:
// g'(x) = exp(x)/Z - exp(x) exp(Z)/(exp(Z) - 1) - 1
template<class Type>
Type VPA_F(Type logF, Type M, Type CAA, Type N_next) {
  Type Z = exp(logF) + M;
  return log(Z) + log(CAA) - log(N_next) - log(exp(Z) - 1) - logF;
}

// Derivative of VPA_F
template<class Type>
Type deriv_VPA_F(Type logF, Type M) {
  Type Z = exp(logF) + M;
  Type ans = exp(logF)/Z;
  ans -= exp(logF) * exp(Z)/(exp(Z) - 1);
  return ans - 1;
}

// Newton solver for log(F)
template<class Type>
Type Newton_VPA_F(Type F_start, Type M, Type CAA, Type N_next, int nloop) {
  Type logF = log(F_start);
  for(int i=0;i<nloop;i++) {
    Type tmp = VPA_F(logF, M, CAA, N_next)/deriv_VPA_F(logF, M);
    logF -= tmp;
  }
  return exp(logF);
}

// Iterative solver for F for A-1 and A (maximum age as a plus-group) in VPA model -
// Solve for x = log(F_y,A-1) such that:
// h(x) = N_y,A exp(-Z_y,A) + N_y,A-1 exp(-Z_y,A-1) - N_y+1,A = 0
// where:
// N_y,A-1 = Z_y,A-1 / F_y,A-1 / (1 - exp(-Z_y,A-1)) * C_y,A-1
// N_y,A = Z_y,A / F_y,A / (1 - exp(-Z_y,A)) * C_y,A where F_y,A = phi * F_y,A-1
template<class Type>
Type VPA_F_plus(Type logF, Type phi, Type M1, Type M2, Type C1, Type C2, Type N_next) {
  Type F1 = exp(logF);
  Type Z1 = F1 + M1;
  Type F2 = phi * F1;
  Type Z2 = F2 + M2;
  Type N1 = Z1/F1/(1 - exp(-Z1)) * C1;
  Type N2 = Z2/F2/(1 - exp(-Z2)) * C2;
  return N1 * exp(-Z1) + N2 * exp(-Z2) - N_next;
}

// Derivative of h(x)
template<class Type>
Type deriv_VPA_F_plus(Type logF, Type phi, Type M1, Type M2, Type C1, Type C2) {
  Type F1 = exp(logF);
  Type Z1 = F1 + M1;
  Type F2 = phi * F1;
  Type Z2 = F2 + M2;

  Type a = Z1/(1 - exp(-Z1));
  Type a_deriv = 1 - exp(-Z1) - Z1 * exp(-Z1);
  a_deriv *= F1;
  a_deriv /= pow(1 - exp(-Z1), 2);

  Type b = C1/F1;
  Type b_deriv = -C1/F1;

  Type cc = Z2/(1 - exp(-Z2));
  Type cc_deriv = 1 - exp(-Z2) - Z2 * exp(-Z2);
  cc_deriv *= F2;
  cc_deriv /= pow(1 - exp(-Z2), 2);

  Type d = C2/F2;
  Type d_deriv = -C2/F2;

  Type ans = a_deriv * b + a * b_deriv;
  ans += cc_deriv * d + cc * d_deriv;
  return ans;
}

// Newton solver for plus-group F
template<class Type>
Type Newton_VPA_F_plus(Type F_start, Type phi, Type M1, Type M2, Type CAA1, Type CAA2, Type N_next, int nloop) {
  Type logF = log(F_start);
  for(int i=0;i<nloop;i++) {
    Type tmp = VPA_F_plus(logF, phi, M1, M2, CAA1, CAA2, N_next);
    tmp /= deriv_VPA_F_plus(logF, phi, M1, M2, CAA1, CAA2);
    logF -= tmp;
  }
  return exp(logF);
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
        Type SP = CppAD::CondExpEq(n, Type(1), -exp(1) * MSY * B_next / K * log(B_next/K),
                                   nterm/(n-1) * MSY * (B_next/K - pow(B_next/K, n))) - F * B_next;
        SP *= dt;
        Catch += F * B_next * dt;
        B_next += SP;
      }

      F *= C_hist/Catch;
      if(i==nitF-1) {
        B(y+1) = B_next;
        Cpred(y) = Catch;
      }
    }
  } else {
    Cpred(y) = C_hist;
    F = CppAD::CondExpLt(1 - C_hist/B(y), Type(0.025), 1 - posfun(1 - C_hist/B(y), Type(0.025), penalty), C_hist/B(y));
    B(y+1) = B(y) + CppAD::CondExpEq(n, Type(1), -exp(1) * MSY * B(y) / K * log(B(y)/K),
      nterm/(n-1) * MSY * (B(y)/K - pow(B(y)/K, n))) - F * B(y);
  }

  return F;
}


