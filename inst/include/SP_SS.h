
//template<class Type>
//Type objective_function<Type>::operator() ()
//{

  DATA_VECTOR(Catch);
  DATA_VECTOR(Index);
  DATA_INTEGER(n_y);

  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(log_B1frac);
  PARAMETER(log_n);
  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_B_dev);

  Type UMSY = 1/(1 + exp(-logit_UMSY));
  Type MSY = exp(log_MSY);
  Type B1frac = exp(log_B1frac);
  Type n = exp(log_n);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);

  Type gamma = pow(n, n/(n-1))/(n-1);
  Type BMSY = MSY/UMSY;
  Type K = BMSY / pow(n, 1/(1-n));
  Type r = MSY * pow(n, n/(n-1)) / K;

  vector<Type> Biomass(n_y+1);
  vector<Type> Ipred(n_y);
  vector<Type> U(n_y);
  vector<Type> relU(n_y);
  vector<Type> relB(n_y+1);

  Biomass(0) = B1frac * K;
  relB(0) = Biomass(0)/BMSY;
  for(int y=0;y<n_y;y++) {
    U(y) = Catch(y)/Biomass(y);
	relU(y) = U(y)/UMSY;
    Type B_test = Biomass(y) + gamma * MSY * (Biomass(y)/K - pow(Biomass(y)/K, n)) - Catch(y);
    Biomass(y+1) = CppAD::CondExpGt(B_test, Type(1e-15), B_test, Type(1e-15));
	  if(y<n_y-1) Biomass(y+1) *= exp(log_B_dev(y) - 0.5 * pow(tau, 2));
    relB(y+1) = Biomass(y+1)/BMSY;
  }

  Type q = calc_q(Index, Biomass);
  for(int y=0;y<n_y;y++) Ipred(y) = q * Biomass(y);

  vector<Type> nll_comp(2);
  nll_comp.setZero();

  for(int y=0;y<n_y;y++) {
    if(Index(y)>0) nll_comp(0) -= dnorm(log(Index(y)), log(Ipred(y)), sigma, true);
	if(y<n_y-1) nll_comp(1) -= dnorm(log_B_dev(y), Type(0), tau, true);
  }

  Type TAC = UMSY * Biomass(n_y);
  Type nll = nll_comp.sum();

  ADREPORT(UMSY);
  ADREPORT(MSY);
  ADREPORT(B1frac);
  ADREPORT(n);
  ADREPORT(q);
  ADREPORT(tau);
  REPORT(sigma);
  REPORT(tau);
  REPORT(gamma);
  REPORT(B1frac);
  REPORT(n);
  REPORT(q);
  REPORT(r);
  REPORT(K);
  REPORT(BMSY);
  REPORT(Ipred);
  REPORT(Biomass);
  REPORT(relB);
  REPORT(U);
  REPORT(relU);
  REPORT(log_B_dev);
  REPORT(TAC);

  return nll;

//}
