
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
  
  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(U_equilibrium);
  
  PARAMETER(vul_50);
  PARAMETER(log_vul_95_offset);
  
  PARAMETER(log_sigma);
  PARAMETER(log_tau);  
  PARAMETER_VECTOR(log_rec_dev);
  
  Type UMSY = 1/(1 + exp(-logit_UMSY));
  Type MSY = exp(log_MSY);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);
  
  // Vulnerability
  vector<Type> vul(max_age);
  Type vul_95 = vul_50 + exp(log_vul_95_offset);
  for(int a=0;a<max_age;a++) {
	Type aa = a;
	vul(a) = 1/(1 + exp(-log(19) * (aa - vul_50)/(vul_95 - vul_50)));
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
    
  Type k1 = deriv_EPR_UMSY/EPR_UMSY;
  Type k2 = deriv_VBPR_UMSY/VBPR_UMSY;  
  Type Arec = (1 - k1 * UMSY + k2 * UMSY)/(EPR_UMSY * (1 + k2 * UMSY));
  Type Brec = (Arec * EPR_UMSY - 1)/(RMSY * EPR_UMSY);
  
  Type CR = Arec * EPR0;
  Type h = CR/(4 + CR);
  
  Type R0 = (Arec * EPR0 - 1)/(Brec * EPR0);
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
  Type R_equilibrium = CppAD::CondExpGt((Arec * EPR_equilibrium - 1)/(Brec * EPR_equilibrium), Type(1e-8), 
                                        (Arec * EPR_equilibrium - 1)/(Brec * EPR_equilibrium), Type(1e-8));
  
  R(0) = R_equilibrium;
  for(int a=0;a<max_age;a++) {
    N(0,a) = R_equilibrium * NPR_equilibrium(a);
	B(0) += N(0,a) * weight(a);
	VB(0) += N(0,a) * weight(a) * vul(a);
	E(0) += N(0,a) * weight(a)* mat(a);
  }
  
  // Loop over all other years  
  for(int y=0;y<n_y;y++) {
	R(y+1) = CppAD::CondExpGt(BH_SR(E(y), h, R0, E0), Type(1e-8), BH_SR(E(y), h, R0, E0), Type(1e-8));
	if(y<n_y-1) R(y+1) *= exp(log_rec_dev(y) - 0.5 * pow(tau, 2));
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
	if(!R_IsNA(asDouble(CAA_hist(y,0)))) nll_comp(1) -= dmultinom(loglike_CAAobs, loglike_CAApred, true);
	
	if(y>0) nll_comp(2) -= dnorm(log_rec_dev(y-1), Type(0), tau, true);
  }
  
  // Very large penalties to likelihood if stock-recruitment parameters (a and b) are negative.
  Type penalty = CppAD::CondExpLe(Arec, Type(0), Type(1e5), Type(0));
  penalty += CppAD::CondExpLe(Brec, Type(0), Type(1e5), Type(0));
  
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
  
  REPORT(vul_50);
  REPORT(vul_95);
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
  
  //REPORT(U_UMSY_final);
  //REPORT(B_BMSY_final);
  //REPORT(B_B0_final);
  //REPORT(E_EMSY_final);
  //REPORT(E_E0_final);
  //REPORT(VB_VBMSY_final);
  //REPORT(VB_VB0_final);
  
  REPORT(log_rec_dev);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(penalty);
  
  return nll;
//}
