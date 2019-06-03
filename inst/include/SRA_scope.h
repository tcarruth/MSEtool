
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace SRA_scope;

  DATA_VECTOR(C_hist);    // Total catch
  DATA_SCALAR(C_eq);
  DATA_MATRIX(I_hist);    // Index by year and fleet
  DATA_INTEGER(nfleet);
  DATA_MATRIX(CAA_hist);  // Catch-at-age proportions
  DATA_VECTOR(CAA_n);     // Annual samples in CAA
  DATA_MATRIX(CAL_hist);  // Catch-at-length proportions
  DATA_VECTOR(CAL_n);     // Annual samples in CAL
  DATA_VECTOR(length_bin);// Vector of length bins
  DATA_VECTOR(mlen);      // Vector of annual mean lengths
  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_MATRIX(M);         // Natural mortality at age
  DATA_MATRIX(len_age);    // Length-at-age
  DATA_SCALAR(CV_LAA);    // CV of length-at-age
  DATA_VECTOR(wt_at_len);    // Weight-at-length
  DATA_MATRIX(mat);       // Maturity-at-age at the beginning of the year
  DATA_STRING(vul_type);  // String indicating whether logistic or dome vul is used
  DATA_IVECTOR(I_type);    // String whether index surveys B, VB, or SSB
  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_STRING(CAA_dist);  // String indicating whether CAA is multinomial or lognormal
  DATA_INTEGER(est_C_eq); //
  DATA_VECTOR(est_early_rec_dev);
  DATA_VECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero

  PARAMETER(log_R0);
  PARAMETER(transformed_h);
  PARAMETER(log_U_equilibrium);
  PARAMETER_VECTOR(vul_par);

  PARAMETER_VECTOR(log_sigma_I);
  PARAMETER(log_sigma_mlen);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_early_rec_dev);
  PARAMETER_VECTOR(log_rec_dev);

  int nlbin = length_bin.size();
  Type bin_width = length_bin(1) - length_bin(0);

  Type R0 = exp(log_R0);
  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else {
    h = exp(transformed_h);
  }
  h += 0.2;
  vector<Type> sigma_I(nfleet);
  for(int ff=0;ff<nfleet;ff++) sigma_I(ff) = exp(log_sigma_I(ff));

  Type sigma_mlen = exp(log_sigma_mlen);
  Type tau = exp(log_tau);
  Type U_equilibrium = exp(log_U_equilibrium);

  Type penalty = 0;
  Type prior = 0.;

  // Vulnerability (length-based)
  vector<Type> vul(nlbin);
  if(vul_type == "logistic") {
    vul = calc_logistic_vul(vul_par, nlbin, length_bin, prior);
  } else {
    vul = calc_dome_vul(vul_par, nlbin, length_bin, prior);
  }

  ////// Equilibrium reference points and per-recruit quantities
  matrix<Type> NPR_unfished(max_age, nlbin);
  matrix<Type> ALK_unfished(max_age, nlbin);
  ALK_unfished = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, n_y-1);
  NPR_unfished = calc_NPR(Type(0), vul, nlbin, M, max_age, ALK_unfished, n_y-1);

  Type EPR0 = sum_EPR(NPR_unfished, wt_at_len, mat, max_age, nlbin, n_y-1);

  Type B0 = R0 * sum_BPR(NPR_unfished, wt_at_len);
  Type N0 = R0 * NPR_unfished.sum();
  Type E0 = R0 * EPR0;
  Type VB0 = R0 * sum_VBPR(NPR_unfished, wt_at_len, vul, max_age, nlbin);

  Type Arec;
  Type Brec;

  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Arec /= EPR0;
    Brec = 5*h - 1;
    Brec /= (1-h) * E0;
  } else {
    Arec = pow(5*h, 1.25);
    Arec /= EPR0;
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= E0;
  }
  Type CR = Arec * EPR0;

  ////// During time series year = 1, 2, ..., n_y
  array<Type> N_full(n_y+1, max_age, nlbin);   // Numbers at year and age
  vector<matrix<Type> > ALK(n_y+1); // sum over lengths = 1
  array<Type> Cat(n_y, max_age, nlbin);
  matrix<Type> N(n_y+1, max_age);
  matrix<Type> CAApred(n_y, max_age);   // Catch (in numbers) at year and age at the mid-point of the season
  matrix<Type> CALpred(n_y, length_bin.size());
  vector<Type> mlen_pred(n_y);
  vector<Type> CN(n_y);             // Catch in numbers
  vector<Type> U(n_y);              // Harvest rate at year
  matrix<Type> Ipred(n_y, nfleet);          // Predicted index at year
  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(max_age-1);
  vector<Type> VB(n_y+1);           // Vulnerable biomass at the midpoint of the year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year

  CN.setZero();
  VB.setZero();
  CAApred.setZero();
  CALpred.setZero();
  mlen_pred.setZero();
  B.setZero();
  E.setZero();
  N.setZero();

  // Equilibrium quantities (leading into first year of model)
  matrix<Type> NPR_equilibrium(max_age, nlbin);
  matrix<Type> ALK_equilibrium(max_age, nlbin);
  ALK_equilibrium = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, 0);
  NPR_equilibrium = calc_NPR(U_equilibrium, vul, nlbin, M, max_age, ALK_equilibrium, 0);

  Type EPR_eq = sum_EPR(NPR_equilibrium, wt_at_len, mat, max_age, nlbin, 0);
  Type R_eq;

  if(SR_type == "BH") {
    R_eq = Arec * EPR_eq - 1;
  } else {
    R_eq = log(Arec * EPR_eq);
  }
  R_eq /= Brec * EPR_eq;

  R(0) = R_eq;
  if(!R_IsNA(asDouble(est_rec_dev(0)))) {
    R(0) *= exp(log_rec_dev(0) - 0.5 * pow(tau, 2));
  }

  ALK(0) = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, 0);
  for(int a=0;a<max_age;a++) {
    if(a == 0) {
      for(int len=0;len<nlbin;len++) N(0,a) += NPR_equilibrium(a,len);
      N(0,a) *= R(0);
    } else {
      R_early(a-1) = R_eq;
      if(!R_IsNA(asDouble(est_early_rec_dev(a-1)))) {
        R_early(a-1) *= exp(log_early_rec_dev(a-1) - 0.5 * pow(tau, 2));
      }
      for(int len=0;len<nlbin;len++) N(0,a) += NPR_equilibrium(a,len);
      N(0,a) *=  R_early(a-1);
    }
    for(int len=0;len<nlbin;len++) {
      N_full(0,a,len) = N(0,a) * ALK(0)(a,len);
      VB(0) += N_full(0,a,len) * wt_at_len(len) * vul(len) * exp(-0.5 * M(0,a));
      B(0) += N_full(0,a,len) * wt_at_len(len);
      E(0) += N_full(0,a,len) * wt_at_len(len) * mat(0,a);
    }
  }
  Type C_eq_pred = calc_C_eq(U_equilibrium, N_full, vul, M, wt_at_len, nlbin, max_age, 0);

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y), h, R0, E0);
    } else {
      R(y+1) = Ricker_SR(E(y), h, R0, E0);
    }

    if(y<n_y-1) {
      if(!R_IsNA(asDouble(est_rec_dev(y)))) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * pow(tau, 2));
    }
    N(y+1,0) = R(y+1);
    U(y) = CppAD::CondExpLt(1 - C_hist(y)/VB(y), Type(0.025),
                            1 - posfun(1 - C_hist(y)/VB(y), Type(0.025), penalty), C_hist(y)/VB(y));

    int x = y+1;
    if(y==n_y-1) x = y;
    ALK(y+1) = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, x);

    for(int a=0;a<max_age;a++) {
      for(int len=0;len<nlbin;len++) {
        Cat(y,a,len) = vul(len) * U(y) * N_full(y,a,len) * exp(-0.5 * M(y,a));
        CAApred(y,a) += Cat(y,a,len);
        CALpred(y,len) += Cat(y,a,len);
        CN(y) += Cat(y,a,len);
        mlen_pred(y) += Cat(y,a,len) * length_bin(len);

        if(a<max_age-1) N(y+1,a+1) += N_full(y,a,len) * exp(-M(y,a)) * (1 - vul(len) * U(y));
        if(a==max_age-1) N(y+1,a) += N_full(y,a,len) * exp(-M(y,a)) * (1 - vul(len) * U(y));

        N_full(y+1,a,len) = N(y+1,a) * ALK(y+1)(a,len);
        VB(y+1) += N_full(y+1,a,len) * wt_at_len(len) * vul(len) * exp(-0.5 * M(y,a));
        B(y+1) += N_full(y+1,a,len) * wt_at_len(len);
        E(y+1) += N_full(y+1,a,len) * wt_at_len(len) * mat(y,a);
      }
    }
    mlen_pred(y) /= CN(y);
  }

  // Calculate nuisance parameters and likelihood
  vector<Type> q(nfleet);
  for(int ff=0;ff<nfleet;ff++) {
    if(I_type(ff) == 0) { // "B"
      q(ff) = calc_q(I_hist, B, ff);
      for(int y=0;y<n_y;y++) Ipred(y,ff) = q(ff) * B(y);
    } else if(I_type(ff) == 1) { // "VB"
      q(ff) = calc_q(I_hist, VB, ff);
      for(int y=0;y<n_y;y++) Ipred(y,ff) = q(ff) * VB(y);
    } else { // "SSB"
      q(ff) = calc_q(I_hist, E, ff);
      for(int y=0;y<n_y;y++) Ipred(y,ff) = q(ff) * E(y);
    }
  }

  // 0 = index, 1 = CAA, 2 = CAL, 3 = mean length, 4 = rec_dev, 5 = eq. catch
  vector<Type> nll_comp(6);
  nll_comp.setZero();
  for(int y=0;y<n_y;y++) {
    for(int ff=0;ff<nfleet;ff++) if(!R_IsNA(asDouble(I_hist(y)))) nll_comp(0) -= dnorm(log(I_hist(y,ff)), log(Ipred(y,ff)), sigma_I(ff), true);
    if(C_hist(y) > 0) {
      vector<Type> loglike_CAAobs(max_age);
      vector<Type> loglike_CAApred(max_age);
      for(int a=0;a<max_age;a++) loglike_CAApred(a) = CAApred(y,a)/CN(y);
      if(!R_IsNA(asDouble(CAA_n(y))) && CAA_n(y) > 0) {
        if(CAA_dist == "multinomial") {
          for(int a=0;a<max_age;a++) loglike_CAAobs(a) = CppAD::CondExpLt(CAA_hist(y,a), Type(1e-8), Type(1e-8), CAA_hist(y,a)) * CAA_n(y);
          nll_comp(1) -= dmultinom(loglike_CAAobs, loglike_CAApred, true);
        } else {
          for(int a=0;a<max_age;a++) loglike_CAAobs(a) = CppAD::CondExpLt(CAA_hist(y,a), Type(1e-8), Type(1e-8), CAA_hist(y,a));
          nll_comp(1) -= dlnorm_comp(loglike_CAAobs, loglike_CAApred);
        }
      }

      vector<Type> loglike_CALobs(nlbin);
      vector<Type> loglike_CALpred(nlbin);
      for(int len=0;len<nlbin;len++) loglike_CALpred(len) = CALpred(y,len)/CN(y);
      if(!R_IsNA(asDouble(CAL_n(y))) && CAL_n(y) > 0) {
        if(CAA_dist == "multinomial") {
          for(int len=0;len<nlbin;len++) loglike_CALobs(len) = CppAD::CondExpLt(CAL_hist(y,len), Type(1e-8), Type(1e-8), CAL_hist(y,len)) * CAL_n(y);
          nll_comp(2) -= dmultinom(loglike_CALobs, loglike_CALpred, true);
        } else {
          for(int len=0;len<nlbin;len++) loglike_CALobs(len) = CppAD::CondExpLt(CAL_hist(y,len), Type(1e-8), Type(1e-8), CAL_hist(y,len));
          nll_comp(2) -= dlnorm_comp(loglike_CALobs, loglike_CALpred);
        }
      }

    }
    if(!R_IsNA(asDouble(mlen(y))) && mlen(y) > 0) nll_comp(3) -= dnorm(mlen(y), mlen_pred(y), sigma_mlen, true);
    if(!R_IsNA(asDouble(est_rec_dev(y)))) nll_comp(4) -= dnorm(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<max_age-1;a++) {
    if(!R_IsNA(asDouble(est_early_rec_dev(a)))) nll_comp(4) -= dnorm(log_early_rec_dev(a), Type(0), tau, true);
  }
  if(est_C_eq > 0) nll_comp(5) -= dnorm(log(C_eq), log(C_eq_pred), Type(0.01), true);

  Type nll = nll_comp.sum() + penalty + prior;

  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(sigma_I);
  ADREPORT(tau);
  ADREPORT(q);

  REPORT(sigma_I);
  REPORT(tau);

  REPORT(NPR_unfished);
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
  REPORT(N_full);
  REPORT(ALK);
  REPORT(Cat);
  REPORT(CN);
  REPORT(CAApred);
  REPORT(CALpred);
  REPORT(mlen_pred);
  REPORT(U);
  REPORT(Ipred);
  REPORT(R);
  REPORT(R_early);
  REPORT(VB);
  REPORT(B);
  REPORT(E);
  REPORT(C_eq_pred);

  REPORT(log_early_rec_dev);
  REPORT(log_rec_dev);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);


  return nll;
//}


