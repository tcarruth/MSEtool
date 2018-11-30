
//template<class Type>
//Type objective_function<Type>::operator() ()
//{

  DATA_VECTOR(C_hist);
  DATA_VECTOR(I_hist);
  DATA_INTEGER(ny);
  DATA_VECTOR(est_B_dev);
  DATA_VECTOR_INDICATOR(keep, I_hist);

  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(log_dep);
  PARAMETER(log_n);
  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_B_dev);

  Type UMSY = invlogit(logit_UMSY);
  Type MSY = exp(log_MSY);
  Type dep = exp(log_dep);
  Type n = exp(log_n);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);

  Type n_term = CppAD::CondExpEq(n, Type(1), exp(1), pow(n, n/(n-1)));
  Type n_term2 = CppAD::CondExpEq(n, Type(1), 1/exp(1), pow(n, 1/(1-n)));

  Type gamma = n_term/(n-1);
  Type BMSY = MSY/UMSY;
  Type K = BMSY / n_term2;
  Type r = MSY * n_term / K;

  vector<Type> B(ny+1);
  vector<Type> SP(ny);
  vector<Type> Ipred(ny);
  vector<Type> U(ny);

  Type penalty = 0;

  B(0) = dep * K;
  if(!R_IsNA(asDouble(est_B_dev(0)))) B(0) *= exp(log_B_dev(0) - 0.5 * pow(tau, 2));
  for(int y=0;y<ny;y++) {
    U(y) = CppAD::CondExpLt(1 - C_hist(y)/B(y), Type(0.025),
      1 - posfun(1 - C_hist(y)/B(y), Type(0.025), penalty), C_hist(y)/B(y));
    SP(y) = gamma * MSY * (B(y)/K - pow(B(y)/K, n));
    B(y+1) = B(y) + SP(y) - U(y) * B(y);
	  if(y<ny-1) {
	    if(!R_IsNA(asDouble(est_B_dev(y+1)))) B(y+1) *= exp(log_B_dev(y+1) - 0.5 * pow(tau, 2));
	  }
  }

  Type q = calc_q(I_hist, B);
  for(int y=0;y<ny;y++) Ipred(y) = q * B(y);

  vector<Type> nll_comp(2);
  nll_comp.setZero();

  for(int y=0;y<ny;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll_comp(0) -= keep(y) * dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
    if(!R_IsNA(asDouble(est_B_dev(y)))) nll_comp(1) -= dnorm(log_B_dev(y), Type(0), tau, true);
  }

  Type nll = nll_comp.sum() + penalty;

  Type U_UMSY_final = U(U.size()-1)/UMSY;
  Type B_BMSY_final = B(B.size()-1)/BMSY;
  Type B_K_final = B(B.size()-1)/K;

  ADREPORT(UMSY);
  ADREPORT(MSY);
  ADREPORT(dep);
  ADREPORT(n);
  ADREPORT(q);
  ADREPORT(r);
  ADREPORT(K);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(U_UMSY_final);
  ADREPORT(B_BMSY_final);
  ADREPORT(B_K_final);
  REPORT(UMSY);
  REPORT(MSY);
  REPORT(dep);
  REPORT(n);
  REPORT(sigma);
  REPORT(tau);
  REPORT(gamma);
  REPORT(r);
  REPORT(K);
  REPORT(BMSY);
  REPORT(Ipred);
  REPORT(B);
  REPORT(SP);
  REPORT(U);
  REPORT(log_B_dev);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(penalty);

  return nll;

//}
