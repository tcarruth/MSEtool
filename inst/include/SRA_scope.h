
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace SRA_scope;

  DATA_MATRIX(C_hist);    // Total catch by year and fleet
  DATA_VECTOR(C_eq);      // Equilibrium catch by fleet

  DATA_MATRIX(I_hist);    // Index by year and survey
  DATA_MATRIX(sigma_I);   // Standard deviation of index by year and survey

  DATA_ARRAY(CAA_hist);   // Catch-at-age proportions by year, age, fleet
  DATA_MATRIX(CAA_n);     // Annual samples in CAA by year and fleet

  DATA_ARRAY(CAL_hist);   // Catch-at-length proportions by year, length_bin, fleet
  DATA_MATRIX(CAL_n);     // Annual samples in CAL by year and fleet

  DATA_VECTOR(length_bin);// Vector of length bins
  DATA_MATRIX(mlen);      // Vector of annual mean lengths by year and fleet

  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_INTEGER(nfleet);   // Number of fleets
  DATA_INTEGER(nsurvey);  // Number of surveys

  DATA_MATRIX(M);         // Natural mortality at age
  DATA_MATRIX(len_age);   // Length-at-age
  DATA_SCALAR(CV_LAA);    // CV of length-at-age
  DATA_VECTOR(wt_at_len); // Weight-at-length
  DATA_MATRIX(mat);       // Maturity-at-age at the beginning of the year

  DATA_IVECTOR(vul_type); // Integer vector indicating whether logistic (0) or dome vul (1) is used
  DATA_IVECTOR(I_type);   // Integer vector indicating the basis of the indices for fleet (1-nfleet) or surveys B (-1) or SSB (-2)

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_MATRIX(LWT_C);     // LIkelihood weights for catch, CAA, CAL, ML, C_eq
  DATA_VECTOR(LWT_Index); // LIkelihood weights for the index

  DATA_VECTOR(est_early_rec_dev);
  DATA_VECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero

  PARAMETER(log_R0);
  PARAMETER(transformed_h);
  PARAMETER_MATRIX(vul_par);            // Matrix of vul_par
  //PARAMETER_VECTOR(log_mean_F);
  PARAMETER_MATRIX(log_F);
  PARAMETER_VECTOR(log_F_equilibrium);  // Equilibrium U by fleet

  PARAMETER_VECTOR(log_sigma_mlen);
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

  Type tau = exp(log_tau);
  vector<Type> sigma_mlen(nfleet);
  vector<Type> F_equilibrium(nfleet);
  F_equilibrium.setZero();
  for(int ff=0;ff<nfleet;ff++) {
    sigma_mlen(ff) = exp(log_sigma_mlen(ff));
    if(C_eq(ff)>0) F_equilibrium(ff) = exp(log_F_equilibrium(ff));
  }

  Type penalty = 0;
  Type prior = 0.;

  // Vulnerability (length-based)
  vector<vector<Type> > vul(nfleet);
  vul = calc_vul(vul_par, vul_type, nlbin, length_bin, prior);

  matrix<Type> F(n_y,nfleet);
  for(int ff=0;ff<nfleet;ff++) {
    //for(int y=0;y<n_y;y++) F(y,ff) = exp(log_mean_F(ff) + log_F_dev(y,ff));
    for(int y=0;y<n_y;y++) F(y,ff) = exp(log_F(y,ff));
  }

  array<Type> Z_total(n_y, max_age, nlbin);
  //Z_total.setZero();

  ////// Equilibrium reference points and per-recruit quantities
  matrix<Type> NPR_unfished(max_age, nlbin);
  matrix<Type> ALK_unfished(max_age, nlbin);
  ALK_unfished = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, n_y-1);
  NPR_unfished = calc_NPR0(nlbin, M, max_age, ALK_unfished, n_y-1);

  Type EPR0 = sum_EPR(NPR_unfished, wt_at_len, mat, max_age, nlbin, n_y-1);

  Type B0 = R0 * sum_BPR(NPR_unfished, wt_at_len);
  Type N0 = R0 * NPR_unfished.sum();
  Type E0 = R0 * EPR0;
  //Type VB0 = R0 * sum_VBPR(NPR_unfished, wt_at_len, vul, max_age, nlbin);

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
  vector<matrix<Type> > ALK(n_y+1);
  array<Type> N_full(n_y+1, max_age, nlbin);   // Numbers at year, age, and length
  matrix<Type> N(n_y+1, max_age);

  array<Type> Cat(n_y, max_age, nlbin);
  array<Type> CAApred(n_y, max_age, nfleet);   // Catch (in numbers) at year and age at the mid-point of the season
  array<Type> CALpred(n_y, nlbin, nfleet);
  matrix<Type> mlen_pred(n_y, nfleet);
  matrix<Type> CN(n_y, nfleet);             // Catch in numbers

  matrix<Type> Cpred(n_y, nfleet);
  matrix<Type> Ipred(n_y, nsurvey);          // Predicted index at year

  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(max_age-1);
  matrix<Type> VB(n_y+1, nfleet);   // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year

  N.setZero();

  Cat.setZero();
  CAApred.setZero();
  CALpred.setZero();
  mlen_pred.setZero();
  CN.setZero();

  Cpred.setZero();
  Ipred.setZero();

  VB.setZero();
  B.setZero();
  E.setZero();

  // Equilibrium quantities (leading into first year of model)
  matrix<Type> NPR_equilibrium(max_age, nlbin);
  matrix<Type> ALK_equilibrium(max_age, nlbin);
  ALK_equilibrium = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, 0);
  NPR_equilibrium = calc_NPR(F_equilibrium, vul, nfleet, nlbin, M, max_age, ALK_equilibrium, 0);

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
      N(0,a) = R(0) * get_N_at_age(NPR_equilibrium, nlbin, a);
    } else {
      R_early(a-1) = R_eq;
      if(!R_IsNA(asDouble(est_early_rec_dev(a-1)))) {
        R_early(a-1) *= exp(log_early_rec_dev(a-1) - 0.5 * pow(tau, 2));
      }
      N(0,a) = R_early(a-1) * get_N_at_age(NPR_equilibrium, nlbin, a);
    }
    for(int len=0;len<nlbin;len++) {
      N_full(0,a,len) = N(0,a) * ALK(0)(a,len);
      for(int ff=0;ff<nfleet;ff++) VB(0,ff) += N_full(0,a,len) * wt_at_len(len) * vul(ff)(len);
      B(0) += N_full(0,a,len) * wt_at_len(len);
      E(0) += N_full(0,a,len) * wt_at_len(len) * mat(0,a);
    }
  }

  vector<Type> C_eq_pred(nfleet);
  C_eq_pred = calc_C_eq(F_equilibrium, N_full, vul, M, wt_at_len, nlbin, nfleet, max_age, 0);

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

    int x = y+1;
    if(y==n_y-1) x = y;
    ALK(y+1) = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, x);

    for(int a=0;a<max_age;a++) {
      for(int len=0;len<nlbin;len++) {
        Z_total(y,a,len) = M(y,a);
        for(int ff=0;ff<nfleet;ff++) Z_total(y,a,len) += vul(ff)(len) * F(y,ff);

        Type mean_N = N_full(y,a,len) * (1 - exp(-Z_total(y,a,len)));
        mean_N /= Z_total(y,a,len);

        for(int ff=0;ff<nfleet;ff++) {
          Type Cat_yalf = vul(ff)(len) * F(y,ff) * mean_N;

          CAApred(y,a,ff) += Cat_yalf;
          CALpred(y,len,ff) += Cat_yalf;
          Cat(y,a,len) += Cat_yalf;
          CN(y,ff) += Cat_yalf;
          mlen_pred(y,ff) += Cat_yalf * length_bin(len);
          Cpred(y,ff) += Cat_yalf * wt_at_len(len);
        }

        if(a<max_age-1) N(y+1,a+1) += N_full(y,a,len) * exp(-Z_total(y,a,len));
        if(a==max_age-1) N(y+1,a) += N_full(y,a,len) * exp(-Z_total(y,a,len));
      }

      for(int len=0;len<nlbin;len++) {
        N_full(y+1,a,len) = N(y+1,a) * ALK(y+1)(a,len);
        for(int ff=0;ff<nfleet;ff++) VB(y+1,ff) += vul(ff)(len) * N_full(y+1,a,len) * wt_at_len(len);
        B(y+1) += N_full(y+1,a,len) * wt_at_len(len);
        E(y+1) += N_full(y+1,a,len) * wt_at_len(len) * mat(x,a);
      }
    }
    for(int ff=0;ff<nfleet;ff++) mlen_pred(y,ff) /= CN(y,ff);
  }

  // Calculate nuisance parameters and likelihood
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    if(I_type(sur) > 0) { // VB.col(sur);
      q(sur) = calc_q(I_hist, VB, sur, I_type(sur), Ipred);
    } else if(I_type(sur) == -1) { // "B"
      q(sur) = calc_q(I_hist, B, sur, Ipred);
    } else {
      q(sur) = calc_q(I_hist, E, sur, Ipred);
    }
  }

  // 0 = index, 1 = CAA, 2 = CAL, 3 = mean length, 4 = rec_dev, 5 = eq. catch
  vector<Type> nll_Catch(nfleet);
  vector<Type> nll_Index(nsurvey);
  vector<Type> nll_CAA(nfleet);
  vector<Type> nll_CAL(nfleet);
  vector<Type> nll_ML(nfleet);
  Type nll_log_rec_dev = 0;
  vector<Type> nll_Ceq(nfleet);

  nll_Catch.setZero();
  nll_Index.setZero();
  nll_CAA.setZero();
  nll_CAL.setZero();
  nll_ML.setZero();
  //nll_Ceq.setZero();

  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      if(!R_IsNA(asDouble(I_hist(y,sur)))) nll_Index(sur) -= dnorm(log(I_hist(y,sur)), log(Ipred(y,sur)), sigma_I(y,sur), true);
    }
    nll_Index(sur) *= LWT_Index(sur);
  }

  for(int ff=0;ff<nfleet;ff++) {
    for(int y=0;y<n_y;y++) {
      if(C_hist(y,ff)>0) {

        nll_Catch(ff) -= dnorm(log(C_hist(y,ff)), log(Cpred(y,ff)), Type(0.01), true);

        if(!R_IsNA(asDouble(CAA_n(y,ff))) && CAA_n(y,ff) > 0) {
          vector<Type> loglike_CAAobs(max_age);
          vector<Type> loglike_CAApred(max_age);
          for(int a=0;a<max_age;a++) loglike_CAApred(a) = CAApred(y,a,ff)/CN(y,ff);

          for(int a=0;a<max_age;a++) loglike_CAAobs(a) = CppAD::CondExpLt(CAA_hist(y,a,ff), Type(1e-8), Type(1e-8), CAA_hist(y,a,ff)) * CAA_n(y,ff);
          nll_CAA(ff) -= dmultinom(loglike_CAAobs, loglike_CAApred, true);
        }

        if(!R_IsNA(asDouble(CAL_n(y,ff))) && CAL_n(y,ff) > 0) {
          vector<Type> loglike_CALobs(nlbin);
          vector<Type> loglike_CALpred(nlbin);
          for(int len=0;len<nlbin;len++) loglike_CALpred(len) = CALpred(y,len,ff)/CN(y,ff);

          for(int len=0;len<nlbin;len++) loglike_CALobs(len) = CppAD::CondExpLt(CAL_hist(y,len,ff), Type(1e-8), Type(1e-8), CAL_hist(y,len,ff)) * CAL_n(y,ff);
          nll_CAL(ff) -= dmultinom(loglike_CALobs, loglike_CALpred, true);
        }

        if(!R_IsNA(asDouble(mlen(y,ff))) && mlen(y,ff) > 0) {
          nll_ML(ff) -= dnorm(mlen(y,ff), mlen_pred(y,ff), sigma_mlen(ff), true);
        }
      }
    }

    nll_Catch(ff) *= LWT_C(ff,0);
    nll_CAA(ff) *= LWT_C(ff,1);
    nll_CAL(ff) *= LWT_C(ff,2);
    nll_ML(ff) *= LWT_C(ff,3);
    if(C_eq(ff) > 0) nll_Ceq(ff) = LWT_C(ff,4) * dnorm(log(C_eq(ff)), log(C_eq_pred(ff)), Type(0.01), true);

  }

  for(int y=0;y<n_y;y++) {
    if(!R_IsNA(asDouble(est_rec_dev(y)))) nll_log_rec_dev -= dnorm(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<max_age-1;a++) {
    if(!R_IsNA(asDouble(est_early_rec_dev(a)))) nll_log_rec_dev -= dnorm(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_Catch.sum() + nll_Index.sum();
  nll += nll_CAA.sum() + nll_CAL.sum() + nll_ML.sum();
  nll += nll_log_rec_dev + nll_Ceq.sum();
  nll += penalty + prior;

  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(tau);
  ADREPORT(q);

  REPORT(M);

  REPORT(log_R0);
  REPORT(transformed_h);
  REPORT(vul_par);
  //REPORT(log_mean_F);
  REPORT(log_F);
  REPORT(log_F_equilibrium);

  REPORT(log_sigma_mlen);
  REPORT(log_tau);
  REPORT(log_early_rec_dev);
  REPORT(log_rec_dev);

  REPORT(R0);
  REPORT(h);
  REPORT(tau);
  REPORT(sigma_mlen);
  REPORT(F_equilibrium);
  REPORT(vul);
  REPORT(F);
  REPORT(Z_total);

  REPORT(NPR_unfished);
  REPORT(ALK_unfished);
  REPORT(EPR0);
  REPORT(B0);
  REPORT(E0);
  REPORT(N0);

  REPORT(Arec);
  REPORT(Brec);
  REPORT(CR);

  REPORT(N_full);
  REPORT(ALK);
  REPORT(Cat);
  REPORT(N);
  REPORT(CAApred);
  REPORT(CALpred);
  REPORT(mlen_pred);
  REPORT(CN);
  REPORT(Cpred);
  REPORT(Ipred);
  REPORT(R);
  REPORT(R_early);
  REPORT(VB);
  REPORT(B);
  REPORT(E);

  REPORT(NPR_equilibrium);
  REPORT(ALK_equilibrium);
  REPORT(EPR_eq);
  REPORT(R_eq);

  REPORT(C_eq_pred);
  REPORT(q);

  REPORT(nll_Catch);
  REPORT(nll_Index);
  REPORT(nll_CAA);
  REPORT(nll_CAL);
  REPORT(nll_ML);
  REPORT(nll_Ceq);

  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  return nll;

//}


