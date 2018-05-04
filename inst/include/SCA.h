//template<class Type>
//Type objective_function<Type>::operator() ()
//{

  DATA_VECTOR(C_hist);    // Total catch
  DATA_VECTOR(I_hist);    // Index
  DATA_MATRIX(CAA_hist);  // Catch-at-age proportions
  DATA_VECTOR(CAA_n);     // Annual samples in CAA
  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_VECTOR(M);         // Natural mortality at age
  DATA_VECTOR(weight);    // Weight-at-age at the beginning of the year
  DATA_VECTOR(mat);       // Maturity-at-age at the beginning of the year
  DATA_STRING(vul_type);  // String indicating whether logistic or dome vul is used
  DATA_VECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero

  PARAMETER(log_meanR);
  PARAMETER(U_equilibrium);
  PARAMETER_VECTOR(vul_par);

  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_rec_dev);

  Type meanR = exp(log_meanR);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);

  // Vulnerability
  vector<Type> vul(max_age);
  if(vul_type == "logistic") {
    vul = calc_logistic_vul(vul_par, max_age);
  } else {
    vul = calc_dome_vul(vul_par, max_age);
  }

  ////// Equilibrium reference points and per-recruit quantities
  vector<Type> NPR_virgin(max_age);
  NPR_virgin = calc_NPR(Type(0), vul, M, max_age);                     // Numbers-per-recruit (NPR) at U = 0
  // Virgin reference points and stock-recruit parameters
  Type EPR0 = sum_EPR(NPR_virgin, weight, mat);                        // Egg-per-recruit at U = 0

  ////// During time series year = 1, 2, ..., n_y
  matrix<Type> N(n_y+1, max_age);   // Numbers at year and age
  matrix<Type> CAApred(n_y, max_age);   // Catch (in numbers) at year and age at the mid-point of the season
  vector<Type> CN(n_y);             // Catch in numbers
  vector<Type> U(n_y);              // Harvest rate at year
  vector<Type> Ipred(n_y);          // Predicted index at year
  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> VB(n_y+1);           // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year

  CN.setZero();
  VB.setZero();
  B.setZero();
  E.setZero();

  // Equilibrium quantities (leading into first year of model)
  vector<Type> NPR_equilibrium(max_age);
  NPR_equilibrium = calc_NPR(U_equilibrium, vul, M, max_age);

  R(0) = meanR;
  if(!R_IsNA(asDouble(est_rec_dev(0)))) {
    R(0) *= exp(log_rec_dev(0) - 0.5 * pow(tau, 2));
  }
  for(int a=0;a<max_age;a++) {
    N(0,a) = R(0) * NPR_equilibrium(a);
    B(0) += N(0,a) * weight(a);
    VB(0) += N(0,a) * weight(a) * vul(a);
    E(0) += N(0,a) * weight(a) * mat(a);
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    R(y+1) = meanR;
    if(y<n_y-1) {
      if(!R_IsNA(asDouble(est_rec_dev(y)))) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * pow(tau, 2));
    }
    N(y+1,0) = R(y+1);

    U(y) = C_hist(y)/VB(y);
    for(int a=0;a<max_age;a++) {
      CAApred(y,a) = vul(a) * U(y) * N(y,a);
      CN(y) += CAApred(y,a);
      if(a<max_age-1) {
        N(y+1,a+1) = CppAD::CondExpGt(N(y,a) * exp(-M(a)) * (1 - vul(a) * U(y)), Type(1e-8),
                                      N(y,a) * exp(-M(a)) * (1 - vul(a) * U(y)), Type(1e-8));
	    }
      if(a==max_age-1) N(y+1,a) += CppAD::CondExpGt(N(y,a) * exp(-M(a)) * (1 - vul(a) * U(y)), Type(1e-8),
	                                                  N(y,a) * exp(-M(a)) * (1 - vul(a) * U(y)), Type(1e-8));
	    B(y+1) += N(y+1,a) * weight(a);
	    VB(y+1) += N(y+1,a) * weight(a) * vul(a);
	    E(y+1) += N(y+1,a) * weight(a) * mat(a);
	  }
  }

  // Calculate nuisance parameters and likelihood
  Type q = calc_q(I_hist, VB);
  for(int y=0;y<n_y;y++) Ipred(y) = q * VB(y);

  vector<Type> nll_comp(3);
  nll_comp.setZero();
  for(int y=0;y<n_y;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll_comp(0) -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);

    vector<Type> loglike_CAAobs(max_age);
	  vector<Type> loglike_CAApred(max_age);
	  for(int a=0;a<max_age;a++) {
	    loglike_CAAobs(a) = (CAA_hist(y,a) + 1e-8) * CAA_n(y);
	    loglike_CAApred(a) = CAApred(y,a)/CN(y);
	  }
	  if(!R_IsNA(asDouble(CAA_n(y)))) nll_comp(1) -= dmultinom(loglike_CAAobs, loglike_CAApred, true);
	  if(!R_IsNA(asDouble(est_rec_dev(y)))) nll_comp(2) -= dnorm(log_rec_dev(y), Type(0), tau, true);
  }

  Type nll = nll_comp.sum();

  ADREPORT(meanR);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(q);

  REPORT(NPR_virgin);
  REPORT(EPR0);
  REPORT(meanR);

  REPORT(vul_par);
  REPORT(vul);

  REPORT(N);
  REPORT(CN);
  REPORT(CAApred);
  REPORT(U);
  REPORT(Ipred);
  REPORT(R);
  REPORT(VB);
  REPORT(B);
  REPORT(E);

  REPORT(log_rec_dev);
  REPORT(nll_comp);
  REPORT(nll);

  return nll;
//}


