
//template<class Type>
//Type objective_function<Type>::operator() ()
//{

  DATA_VECTOR(C_hist);
  DATA_VECTOR(I_hist);
  DATA_INTEGER(ny);
  DATA_INTEGER(nstep);
  DATA_SCALAR(dt);
  DATA_INTEGER(nitF);
  DATA_VECTOR(r_prior);

  PARAMETER(log_FMSY);
  PARAMETER(log_MSY);
  PARAMETER(log_dep);
  PARAMETER(log_n);

  Type FMSY = exp(log_FMSY);
  Type MSY = exp(log_MSY);
  Type dep = exp(log_dep);
  Type n = exp(log_n);

  Type euler = exp(Type(1.0));

  Type n_term = CppAD::CondExpEq(n, Type(1), euler, pow(n, n/(n-1)));
  Type n_term2 = CppAD::CondExpEq(n, Type(1), Type(1/euler), pow(n, 1/(1-n)));

  Type BMSY = MSY/FMSY;
  Type K = BMSY / n_term2;
  Type r = MSY * n_term / K; // r = FMSY * pow(n, 1/(1-n))

  vector<Type> B(ny+1);
  vector<Type> SP(ny);
  vector<Type> Cpred(ny);
  vector<Type> Ipred(ny);
  vector<Type> F(ny);

  Type nll = 0;
  Type penalty = 0;
  Type prior = 0;

  if(r_prior(0) > 0) prior -= dnorm(r, r_prior(0), r_prior(1), true) + log_FMSY; // r prior with log-Jacobian transformation, exact with fixed n

  Cpred.setZero();

  B(0) = dep * K;
  for(int y=0;y<ny;y++) {
    F(y) = SP_F(C_hist(y)/(C_hist(y) + B(y)), C_hist(y), MSY, K, n, n_term,
      CppAD::CondExpLe(C_hist(y), Type(1e-8), Type(1), dt), nstep, nitF, Cpred, B, y, penalty);
    SP(y) = B(y+1) - B(y) + Cpred(y);
  }

  Type q = calc_q(I_hist, B);
  for(int y=0;y<ny;y++) Ipred(y) = q * B(y);

  Type sigma = calc_sigma(I_hist, Ipred);
  for(int y=0;y<ny;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
  }

  nll += penalty + prior;

  Type F_FMSY_final = F(F.size()-1)/FMSY;
  Type B_BMSY_final = B(B.size()-1)/BMSY;
  Type B_K_final = B(B.size()-1)/K;

  ADREPORT(FMSY);
  ADREPORT(MSY);
  ADREPORT(dep);
  ADREPORT(n);
  ADREPORT(q);
  ADREPORT(r);
  ADREPORT(K);
  ADREPORT(sigma);
  ADREPORT(F_FMSY_final);
  ADREPORT(B_BMSY_final);
  ADREPORT(B_K_final);
  REPORT(FMSY);
  REPORT(MSY);
  REPORT(dep);
  REPORT(n);
  REPORT(sigma);
  REPORT(r);
  REPORT(K);
  REPORT(BMSY);
  REPORT(Cpred);
  REPORT(Ipred);
  REPORT(B);
  REPORT(SP);
  REPORT(F);
  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  return nll;

//}
