
//template<class Type>
//Type objective_function<Type>::operator() ()
//{

  DATA_VECTOR(C_hist);
  DATA_VECTOR(I_hist);
  DATA_INTEGER(ny);

  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(log_dep);
  PARAMETER(log_n);

  Type UMSY = ilogit(logit_UMSY);
  Type MSY = exp(log_MSY);
  Type dep = exp(log_dep);
  Type n = exp(log_n);

  Type gamma = pow(n, n/(n-1))/(n-1);
  Type BMSY = MSY/UMSY;
  Type K = BMSY / pow(n, 1/(1-n));
  Type r = MSY * pow(n, n/(n-1)) / K;

  vector<Type> B(ny+1);
  vector<Type> SP(ny);
  vector<Type> Ipred(ny);
  vector<Type> U(ny);

  Type nll = 0.;
  Type penalty = 0.;

  B(0) = dep * K;
  for(int y=0;y<ny;y++) {
    U(y) = CppAD::CondExpLt(1 - C_hist(y)/B(y), Type(0.025),
      1 - posfun(1 - C_hist(y)/B(y), Type(0.025), penalty), C_hist(y)/B(y));
    SP(y) = gamma * MSY * (B(y)/K - pow(B(y)/K, n));
    B(y+1) = B(y) + SP(y) - U(y) * B(y);
  }

  Type q = calc_q(I_hist, B);
  for(int y=0;y<ny;y++) Ipred(y) = q * B(y);

  Type sigma = calc_sigma(I_hist, Ipred);
  for(int y=0;y<ny;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
  }

  nll += penalty;

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
  ADREPORT(U_UMSY_final);
  ADREPORT(B_BMSY_final);
  ADREPORT(B_K_final);
  REPORT(UMSY);
  REPORT(MSY);
  REPORT(dep);
  REPORT(n);
  REPORT(sigma);
  REPORT(gamma);
  REPORT(r);
  REPORT(K);
  REPORT(BMSY);
  REPORT(Ipred);
  REPORT(B);
  REPORT(SP);
  REPORT(U);
  REPORT(nll);
  REPORT(penalty);

  return nll;

//}
