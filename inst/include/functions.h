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
vector<Type> calc_NPR(Type U, vector<Type> vul, vector<Type> M, int max_age) {
  vector<Type> NPR(max_age);
  NPR(0) = 1.;
  for(int a=1;a<max_age;a++) {
    NPR(a) = NPR(a-1) * exp(-M(a-1)) * (1 - vul(a-1) * U);
  }
  NPR(max_age-1) /= 1 - exp(-M(max_age-1)) * (1 - vul(max_age-1) * U); // Plus-group
  return NPR;
}

template<class Type>
vector<Type> calc_deriv_NPR(Type U, vector<Type> NPR, vector<Type> vul, vector<Type> M, int max_age) {
  vector<Type> deriv_NPR(max_age);
  deriv_NPR(0) = 0.;
  for(int a=1;a<max_age;a++) {
    deriv_NPR(a) = deriv_NPR(a-1) * exp(-M(a-1)) * (1 - vul(a-1) * U) - NPR(a-1) * exp(-M(a-1)) * vul(a-1);
  }
  // Additional calculation for plus-group
  deriv_NPR(max_age-1) *= 1 - exp(-M(max_age-1)) * (1 - vul(max_age-1) * U);
  deriv_NPR(max_age-1) -= NPR(max_age-2) * exp(-M(max_age-2)) * (1 - vul(max_age-2) * U) * exp(-M(max_age-1)) * vul(max_age-1);
  deriv_NPR(max_age-1) /= pow(1 - exp(-M(max_age-1)) * (1 - vul(max_age-1) * U), 2);
  return deriv_NPR;
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
  Type a_full = invlogit(vul_par(0)) * maxage * 0.75;
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

  return vul;
}

template<class Type>
Type dlnorm_comp(vector<Type> obs, vector<Type> pred) {
  Type log_lik = 0.;
  for(int a=0;a<obs.size();a++) log_lik += dnorm(log(obs(a)), log(pred(a)), pow(0.01/pred(a), 0.5), true);
  return log_lik;
}
