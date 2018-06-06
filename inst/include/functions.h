//invlogit
template<class Type>
Type ilogit(Type x) {
  return pow(1 + exp(-x), -1);
}

//posfun from ADMB
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 100 * pow(x - eps, 2));
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
vector<Type> calc_logistic_vul(vector<Type> vul_par, int max_age) {
  vector<Type> vul(max_age);
  Type vul_50 = vul_par(0);
  Type vul_95 = vul_50 + exp(vul_par(1));

	for(int a=0;a<max_age;a++) {
	  Type aa = a;
	  vul(a) = 1/(1 + exp(-log(19) * (aa - vul_50)/(vul_95 - vul_50)));
  }
	return vul;
}

template<class Type>
vector<Type> calc_dome_vul(vector<Type> vul_par, int max_age) {
  vector<Type> vul(max_age);
  Type vul_sd_asc = exp(vul_par(0));
  Type vul_mu_asc = vul_par(1);
  Type vul_mu_des = vul_mu_asc + exp(vul_par(2));
  Type vul_sd_des = exp(vul_par(3));

  Type denom_asc = dnorm(vul_mu_asc, vul_mu_asc, vul_sd_asc, false);
  Type denom_des = dnorm(vul_mu_des, vul_mu_des, vul_sd_des, false);

  for(int a=0;a<max_age;a++) {
    Type aa = a;
    Type vul_asc = dnorm(aa, vul_mu_asc, vul_sd_asc, false);
    vul_asc /= denom_asc;
    Type vul_des = dnorm(aa, vul_mu_des, vul_sd_des, false);
    vul_des /= denom_des;

    vul(a) = CppAD::CondExpLe(aa, vul_mu_asc, vul_asc, CppAD::CondExpLe(aa, vul_mu_des, Type(1), vul_des));
  }
  return vul;
}
