
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
  DATA_STRING(I_type);    // String whether index surveys B, VB, or SSB
  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_VECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero

  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(U_equilibrium);

  PARAMETER_VECTOR(vul_par);

  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_rec_dev);

  Type UMSY = ilogit(logit_UMSY);
  Type MSY = exp(log_MSY);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);

  Type penalty = 0.;

  // Vulnerability
  vector<Type> vul(max_age);
  if(vul_type == "logistic") {
    vul = calc_logistic_vul(vul_par, max_age);
  } else {
    vul = calc_dome_vul(vul_par, max_age);
  }

  ////// Equilibrium reference points and per-recruit quantities
  vector<Type> NPR_virgin(max_age);
  vector<Type> NPR_UMSY(max_age);
  vector<Type> deriv_NPR_UMSY(max_age);

  NPR_virgin = calc_NPR(Type(0), vul, M, max_age);                     // Numbers-per-recruit (NPR) at U = 0
  NPR_UMSY = calc_NPR(UMSY, vul, M, max_age);                          // Numbers-per-recruit at U = UMSY
  deriv_NPR_UMSY = calc_deriv_NPR(UMSY, NPR_UMSY, vul, M, max_age);    // Derivative of NPR wrt to U at UMSY
  Type EPR_UMSY = sum_EPR(NPR_UMSY, weight, mat);                      // Egg-per-recruit at UMSY
  Type VBPR_UMSY = sum_VBPR(NPR_UMSY, weight, vul);                    // Vulnerable biomass-per-recruit at UMSY
  Type deriv_EPR_UMSY = sum_EPR(deriv_NPR_UMSY, weight, mat);          // Derivative of egg-per-recruit at UMSY
  Type deriv_VBPR_UMSY = sum_VBPR(deriv_NPR_UMSY, weight, vul);        // Derivative of vulnerable biomass-per-recruit at UMSY

  // MSY reference points
  Type VBMSY = MSY/UMSY;
  Type RMSY = VBMSY/VBPR_UMSY;
  Type BMSY = RMSY * sum_BPR(NPR_UMSY, weight);
  Type EMSY = RMSY * EPR_UMSY;

  // Virgin reference points and stock-recruit parameters
  Type EPR0 = sum_EPR(NPR_virgin, weight, mat);

  Type Arec;
  Type Brec;
  Type CR;
  Type h;
  Type R0;

  if(SR_type == "BH") {
    Type k1 = deriv_EPR_UMSY/EPR_UMSY;
    Type k2 = deriv_VBPR_UMSY/VBPR_UMSY;

    Arec = (1 - k1 * UMSY + k2 * UMSY)/(EPR_UMSY * (1 + k2 * UMSY));
    Brec = (Arec * EPR_UMSY - 1)/(RMSY * EPR_UMSY);
    CR = Arec * EPR0;
    h = CR/(4 + CR);
    R0 = (Arec * EPR0 - 1)/(Brec * EPR0);
  } else {
    Type num = -1 * UMSY * VBPR_UMSY * deriv_EPR_UMSY;
    Type denom = EPR_UMSY;
    denom *= UMSY * deriv_VBPR_UMSY + VBPR_UMSY;
    denom += num;
    Type log_aphi = num/denom;

    Arec = exp(log_aphi)/EPR_UMSY;
    Brec = log_aphi/RMSY/EPR_UMSY;
    CR = Arec * EPR0;
    h = exp(0.8 * Brec * EPR0);
    h *= 0.2;
    R0 = log(Arec * EPR0);
    R0 /= Brec;
    R0 /= EPR0;
  }
  Type B0 = R0 * sum_BPR(NPR_virgin, weight);
  Type N0 = R0 * NPR_virgin.sum();
  Type E0 = R0 * EPR0;
  Type VB0 = R0 * sum_VBPR(NPR_virgin, weight, vul);

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
  Type EPR_equilibrium = sum_EPR(NPR_equilibrium, weight, mat);
  Type R_equilibrium;
  if(SR_type == "BH") {
    R_equilibrium = CppAD::CondExpGt((Arec * EPR_equilibrium - 1)/(Brec * EPR_equilibrium), Type(1e-8),
                                     (Arec * EPR_equilibrium - 1)/(Brec * EPR_equilibrium), Type(1e-8));
  } else {
    R_equilibrium = CppAD::CondExpGt(log(Arec * EPR_equilibrium)/(Brec * EPR_equilibrium), Type(1e-8),
                                     log(Arec * EPR_equilibrium)/(Brec * EPR_equilibrium), Type(1e-8));
  }

  R(0) = R_equilibrium;
  for(int a=0;a<max_age;a++) {
    N(0,a) = R_equilibrium * NPR_equilibrium(a);
    B(0) += N(0,a) * weight(a);
    VB(0) += N(0,a) * weight(a) * vul(a);
    E(0) += N(0,a) * weight(a)* mat(a);
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    if(SR_type == "BH") {
      R(y+1) = CppAD::CondExpGt(BH_SR(E(y), h, R0, E0), Type(1e-8), BH_SR(E(y), h, R0, E0), Type(1e-8));
    } else {
      R(y+1) = CppAD::CondExpGt(Ricker_SR(E(y), h, R0, E0), Type(1e-8), Ricker_SR(E(y), h, R0, E0), Type(1e-8));
    }

    if(y<n_y-1) {
      if(!R_IsNA(asDouble(est_rec_dev(y)))) R(y+1) *= exp(log_rec_dev(y) - 0.5 * pow(tau, 2));
    }
    N(y+1,0) = R(y+1);

    U(y) = CppAD::CondExpLt(1 - C_hist(y)/VB(y), Type(0.025),
      1 - posfun(1 - C_hist(y)/VB(y), Type(0.025), penalty), C_hist(y)/VB(y));
    for(int a=0;a<max_age;a++) {
      CAApred(y,a) = vul(a) * U(y) * N(y,a);
      CN(y) += CAApred(y,a);
      if(a<max_age-1) N(y+1,a+1) = N(y,a) * exp(-M(a)) * (1 - vul(a) * U(y));
      if(a==max_age-1) N(y+1,a) += N(y,a) * exp(-M(a)) * (1 - vul(a) * U(y));
	    B(y+1) += N(y+1,a) * weight(a);
	    VB(y+1) += N(y+1,a) * weight(a) * vul(a);
	    E(y+1) += N(y+1,a) * weight(a) * mat(a);
	  }
  }

  // Calculate nuisance parameters and likelihood
  Type q;
  if(I_type == "B") {
    q = calc_q(I_hist, B);
    for(int y=0;y<n_y;y++) Ipred(y) = q * B(y);
  } else if(I_type == "VB") {
    q = calc_q(I_hist, VB);
    for(int y=0;y<n_y;y++) Ipred(y) = q * VB(y);
  } else {
    q = calc_q(I_hist, E);
    for(int y=0;y<n_y;y++) Ipred(y) = q * E(y);
  }

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
	  if(y>0) {
	    if(!R_IsNA(asDouble(est_rec_dev(y-1)))) nll_comp(2) -= dnorm(log_rec_dev(y-1), Type(0), tau, true);
	  }
  }

  // Very large penalties to likelihood if stock-recruitment parameters (a and b) are negative.
  // For both B-H and Ricker, Brec > 0 if Arec * EPR_UMSY - 1 > 0.
  // FOr Ricker, Arec is always > 0. Need to find conditions for Arec with B-H S-R relationship.
  penalty += CppAD::CondExpGt(Arec * EPR_UMSY - 1, Type(0), Type(0), Type(UMSY * 1e3));

  Type nll = nll_comp.sum() + penalty;

  Type U_UMSY_final = U(U.size()-1)/UMSY;
  Type B_BMSY_final = B(B.size()-1)/BMSY;
  Type B_B0_final = B(B.size()-1)/B0;
  Type E_EMSY_final = E(E.size()-1)/EMSY;
  Type E_E0_final = E(E.size()-1)/E0;
  Type VB_VBMSY_final = VB(VB.size()-1)/VBMSY;
  Type VB_VB0_final = VB(VB.size()-1)/VB0;

  ADREPORT(UMSY);
  ADREPORT(MSY);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(q);
  ADREPORT(U_UMSY_final);
  ADREPORT(B_BMSY_final);
  ADREPORT(B_B0_final);
  ADREPORT(E_EMSY_final);
  ADREPORT(E_E0_final);
  ADREPORT(VB_VBMSY_final);
  ADREPORT(VB_VB0_final);

  REPORT(NPR_virgin);
  REPORT(NPR_UMSY);
  REPORT(deriv_NPR_UMSY);
  REPORT(EPR_UMSY);
  REPORT(VBPR_UMSY);
  REPORT(deriv_EPR_UMSY);
  REPORT(deriv_VBPR_UMSY);
  REPORT(UMSY);
  REPORT(MSY);
  REPORT(VBMSY);
  REPORT(RMSY);
  REPORT(BMSY);
  REPORT(EMSY);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(EPR0);
  REPORT(CR);
  REPORT(h);
  REPORT(R0);
  REPORT(B0);
  REPORT(N0);
  REPORT(E0);
  REPORT(VB0);

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
  REPORT(penalty);

  return nll;
//}
