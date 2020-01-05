
#ifndef SRA_scope_hpp
#define SRA_scope_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type SRA_scope(objective_function<Type> *obj) {

  using namespace ns_SRA_scope;

  DATA_MATRIX(C_hist);    // Total catch by year and fleet
  DATA_VECTOR(C_eq);      // Equilibrium catch by fleet

  DATA_MATRIX(E_hist);    // Effort by year and fleet
  DATA_VECTOR(E_eq);      // Equilibrium effort by fleet

  DATA_STRING(condition); // Indicates whether the model will condition on effort or catch
  DATA_INTEGER(nll_C);    // Indicates whether there is a likelihood for the catch if condition = "effort"

  DATA_MATRIX(I_hist);    // Index by year and survey
  DATA_MATRIX(sigma_I);   // Standard deviation of index by year and survey

  DATA_ARRAY(CAA_hist);   // Catch-at-age re-weighted by year, age, fleet
  DATA_MATRIX(CAA_n);     // Annual samples in CAA by year and fleet

  DATA_ARRAY(CAL_hist);   // Catch-at-length re-weighted by year, length_bin, fleet
  DATA_MATRIX(CAL_n);     // Annual samples in CAL by year and fleet

  DATA_ARRAY(s_CAA_hist);   // Catch-at-age re-weighted by year, age, survey
  DATA_MATRIX(s_CAA_n);     // Annual samples in CAA by year and survey

  DATA_ARRAY(s_CAL_hist);   // Catch-at-length re-weighted by year, length_bin, survey
  DATA_MATRIX(s_CAL_n);     // Annual samples in CAL by year and survey

  DATA_VECTOR(length_bin);// Vector of length bins
  DATA_MATRIX(mlen);      // Vector of annual mean lengths by year and fleet

  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_INTEGER(nfleet);   // Number of fleets
  DATA_INTEGER(nsurvey);  // Number of surveys

  DATA_MATRIX(M);         // Natural mortality at age
  DATA_MATRIX(len_age);   // Length-at-age
  DATA_SCALAR(Linf);      // Linf
  DATA_SCALAR(CV_LAA);    // CV of length-at-age
  DATA_MATRIX(wt);        // Weight-at-age
  DATA_MATRIX(mat);       // Maturity-at-age at the beginning of the year

  DATA_IVECTOR(vul_type); // Integer vector indicating whether logistic (-1) or dome vul (0), or age-specific (1 - maxage) is used
  DATA_IVECTOR(s_vul_type); // Same but for surveys
  DATA_IVECTOR(I_type);   // Integer vector indicating the basis of the indices for fleet (1-nfleet) or surveys B (-1) or SSB (-2) or estimated (0)
  DATA_IVECTOR(abs_I);    // Boolean, whether index is an absolute (fix q = 1) or relative terms (estimate q)
  DATA_IVECTOR(I_basis);  // Boolean, whether index is biomass based (= 1) or abundance-based (0)

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_MATRIX(LWT_C);     // LIkelihood weights for catch, CAA, CAL, ML, C_eq
  DATA_VECTOR(LWT_Index); // LIkelihood weights for the index
  DATA_STRING(comp_like);

  DATA_SCALAR(max_F);     // Maximum F in the model
  DATA_SCALAR(rescale);   // Catch rescaler, needed in case q = 1
  DATA_INTEGER(ageM);     // Age of maturity used for averaging E0 and EPR0

  DATA_IVECTOR(est_early_rec_dev);
  DATA_IVECTOR(est_rec_dev); // Indicator whether to estimate rec_dev
  DATA_IVECTOR(yindF);

  PARAMETER(log_R0);
  PARAMETER(transformed_h);
  PARAMETER_MATRIX(vul_par);            // Matrix of vul_par
  PARAMETER_MATRIX(s_vul_par);
  PARAMETER_VECTOR(log_q_effort);
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

  // Vulnerability (length-based) and F parameters
  Type penalty = 0;
  Type prior = 0.;
  vector<Type> LFS(nfleet);
  vector<Type> L5(nfleet);
  vector<Type> Vmaxlen(nfleet);
  array<Type> vul = calc_vul(vul_par, vul_type, len_age, LFS, L5, Vmaxlen, Linf);

  vector<Type> q_effort(nfleet);
  vector<Type> sigma_mlen(nfleet);
  vector<Type> F_equilibrium(nfleet);
  F_equilibrium.setZero();

  matrix<Type> F(n_y,nfleet);
  matrix<Type> Z = M;

  for(int ff=0;ff<nfleet;ff++) {
    sigma_mlen(ff) = exp(log_sigma_mlen(ff));
    q_effort(ff) = exp(log_q_effort(ff));
    if(condition == "catch" && C_eq(ff)>0) F_equilibrium(ff) = exp(log_F_equilibrium(ff));
    if(condition == "effort" && E_eq(ff)>0) F_equilibrium(ff) = q_effort(ff) * E_eq(ff);
    if(condition == "catch") {
      Type tmp = max_F - exp(log_F(yindF(ff),ff));
      F(yindF(ff),ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), exp(log_F(yindF(ff),ff)));

      for(int y=0;y<n_y;y++) {
        if(y != yindF(ff)) {
          Type Ftmp = F(yindF(ff),ff) * exp(log_F(y,ff));
          Type tmp2 = max_F - Ftmp;
          F(y,ff) = CppAD::CondExpLt(tmp2, Type(0), max_F - posfun(tmp2, Type(0), penalty), Ftmp);
        }
      }

    } else {
      for(int y=0;y<n_y;y++) {
        Type tmp = max_F - q_effort(ff) * E_hist(y,ff);
        F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty),
          q_effort(ff) * E_hist(y,ff));
      }
    }
  }


  ////// Equilibrium reference points and per-recruit quantities - calculate annually
  vector<vector<Type> > NPR_unfished(n_y);
  vector<Type> EPR0(n_y);
  vector<Type> E0(n_y);
  vector<Type> B0(n_y);
  vector<Type> N0(n_y);

  Type E0_SR = 0;
  for(int y=0;y<n_y;y++) {
    NPR_unfished(y) = calc_NPR0(M, max_age, y);

    EPR0(y) = sum_EPR(NPR_unfished(y), wt, mat, max_age, y);
    E0(y) = R0 * EPR0(y);
    B0(y) = R0 * sum_BPR(NPR_unfished(y), wt, max_age, y);
    N0(y) = R0 * NPR_unfished(y).sum();

    if(y < ageM) E0_SR += E0(y);
  }
  E0_SR /= Type(ageM);
  Type EPR0_SR = E0_SR/R0;

  Type Arec, Brec;
  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Brec = 5*h - 1;
    Brec /= (1-h);
  } else {
    Arec = pow(5*h, 1.25);
    Brec = 1.25;
    Brec *= log(5*h);
  }
  Arec /= EPR0_SR;
  Brec /= E0_SR;

  ////// During time series year = 1, 2, ..., n_y
  vector<matrix<Type> > ALK(n_y);
  matrix<Type> N(n_y+1, max_age);

  vector<Type> C_eq_pred(nfleet);
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

  C_eq_pred.setZero();
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
  vector<Type> NPR_equilibrium = calc_NPR(F_equilibrium, vul, nfleet, M, max_age, 0);
  Type EPR_eq = sum_EPR(NPR_equilibrium, wt, mat, max_age, 0);
  Type R_eq;

  if(SR_type == "BH") {
    R_eq = Arec * EPR_eq - 1;
  } else {
    R_eq = log(Arec * EPR_eq);
  }
  R_eq /= Brec * EPR_eq;

  R(0) = R_eq;
  if(est_rec_dev(0)) R(0) *= exp(log_rec_dev(0) - 0.5 * tau * tau);

  for(int a=0;a<max_age;a++) {
    if(a == 0) {
      N(0,a) = R(0) * NPR_equilibrium(a);
    } else {
      R_early(a-1) = R_eq;
      if(est_early_rec_dev(a-1)) R_early(a-1) *= exp(log_early_rec_dev(a-1) - 0.5 * tau * tau);
      N(0,a) = R_early(a-1) * NPR_equilibrium(a);
    }

    B(0) += N(0,a) * wt(0,a);
    E(0) += N(0,a) * wt(0,a) * mat(0,a);

    Type Z_eq = M(0,a);
    for(int ff=0;ff<nfleet;ff++) Z_eq += vul(0,a,ff) * F_equilibrium(ff);
    Type mean_N_eq = N(0,a) * (1 - exp(-Z_eq)) / Z_eq;
    for(int ff=0;ff<nfleet;ff++) {
      C_eq_pred(ff) += vul(0,a,ff) * F_equilibrium(ff) * mean_N_eq * wt(0,a);
      VB(0,ff) += N(0,a) * wt(0,a) * vul(0,a,ff);
    }
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y), h, R0, E0_SR);
    } else {
      R(y+1) = Ricker_SR(E(y), h, R0, E0_SR);
    }

    if(y<n_y-1) {
      if(est_rec_dev(y+1)) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * tau * tau);
    }
    N(y+1,0) = R(y+1);
    ALK(y) = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, y);

    for(int a=0;a<max_age;a++) {
      for(int ff=0;ff<nfleet;ff++) Z(y,a) += vul(y,a,ff) * F(y,ff);
      Type mean_N = N(y,a) * (1 - exp(-Z(y,a))) / Z(y,a);

      if(a<max_age-1) N(y+1,a+1) = N(y,a) * exp(-Z(y,a));
      if(a==max_age-1) N(y+1,a) += N(y,a) * exp(-Z(y,a));

      for(int ff=0;ff<nfleet;ff++) {
        CAApred(y,a,ff) = vul(y,a,ff) * F(y,ff) * mean_N;
        CN(y,ff) += CAApred(y,a,ff);
        Cpred(y,ff) += CAApred(y,a,ff) * wt(y,a);

        for(int len=0;len<nlbin;len++) {
          CALpred(y,len,ff) += CAApred(y,a,ff) * ALK(y)(a,len);
          mlen_pred(y,ff) += CAApred(y,a,ff) * ALK(y)(a,len) * length_bin(len);
        }

        VB(y+1,ff) += vul(y+1,a,ff) * N(y+1,a) * wt(y+1,a);
      }

      B(y+1) += N(y+1,a) * wt(y+1,a);
      E(y+1) += N(y+1,a) * wt(y+1,a) * mat(y+1,a);
    }
    for(int ff=0;ff<nfleet;ff++) mlen_pred(y,ff) /= CN(y,ff);
  }

  // Calculate nuisance parameters and likelihood
  // Survey selectivity
  vector<Type> s_LFS(nsurvey);
  vector<Type> s_L5(nsurvey);
  vector<Type> s_Vmaxlen(nsurvey);

  array<Type> s_CAApred(n_y, max_age, nsurvey);
  array<Type> s_CALpred(n_y, nlbin, nsurvey);
  matrix<Type> s_CN(n_y, nsurvey);
  matrix<Type> B_sur(n_y, nsurvey); // Biomass vulnerable to the survey

  s_CAApred.setZero();
  s_CALpred.setZero();
  s_CN.setZero();
  B_sur.setZero();

  array<Type> s_vul = calc_vul_sur(s_vul_par, s_vul_type, len_age, s_LFS, s_L5, s_Vmaxlen, Linf, mat, I_type, vul);
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      for(int a=0;a<max_age;a++) {
        s_CAApred(y,a,sur) = s_vul(y,a,sur) * N(y,a);
        s_CN(y,sur) += s_CAApred(y,a,sur);

        if(I_basis(sur)) {
          B_sur(y,sur) += s_CAApred(y,a,sur) * wt(y,a);
        }
        if(!R_IsNA(asDouble(s_CAL_n(y,sur))) && s_CAL_n(y,sur) > 0) {
          for(int len=0;len<nlbin;len++) s_CALpred(y,len,sur) += s_CAApred(y,a,sur) * ALK(y)(a,len);
        }
      }
    }
    if(!I_basis(sur)) B_sur.col(sur) = s_CN.col(sur);
    q(sur) = calc_q(I_hist, B_sur, sur, sur, Ipred, abs_I, rescale);
  }

  vector<Type> nll_Catch(nfleet);
  vector<Type> nll_Index(nsurvey);
  vector<Type> nll_s_CAA(nsurvey);
  vector<Type> nll_s_CAL(nsurvey);
  vector<Type> nll_CAA(nfleet);
  vector<Type> nll_CAL(nfleet);
  vector<Type> nll_ML(nfleet);
  Type nll_log_rec_dev = 0;
  vector<Type> nll_Ceq(nfleet);

  nll_Catch.setZero();
  nll_Index.setZero();
  nll_CAA.setZero();
  nll_CAL.setZero();
  nll_s_CAA.setZero();
  nll_s_CAL.setZero();
  nll_ML.setZero();
  nll_Ceq.setZero();

  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      if(!R_IsNA(asDouble(I_hist(y,sur)))) nll_Index(sur) -= dnorm(log(I_hist(y,sur)), log(Ipred(y,sur)), sigma_I(y,sur), true);

      if(!R_IsNA(asDouble(s_CAA_n(y,sur))) && s_CAA_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_s_CAA(sur) -= comp_multinom(s_CAA_hist, s_CAApred, s_CN, s_CAA_n, y, max_age, sur);
        } else {
          nll_s_CAA(sur) -= comp_lognorm(s_CAA_hist, s_CAApred, s_CN, s_CAA_n, y, max_age, sur);
        }
      }

      if(!R_IsNA(asDouble(s_CAL_n(y,sur))) && s_CAL_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_s_CAL(sur) -= comp_multinom(s_CAL_hist, s_CALpred, s_CN, s_CAL_n, y, nlbin, sur);
        } else {
          nll_s_CAL(sur) = comp_lognorm(s_CAL_hist, s_CALpred, s_CN, s_CAL_n, y, nlbin, sur);
        }
      }
    }
    nll_Index(sur) *= LWT_Index(sur,0);
    nll_s_CAA(sur) *= LWT_Index(sur,1);
    nll_s_CAL(sur) *= LWT_Index(sur,2);
  }

  for(int ff=0;ff<nfleet;ff++) {
    for(int y=0;y<n_y;y++) {
      if(C_hist(y,ff)>0 || E_hist(y,ff)>0) {
        if(!R_IsNA(asDouble(CAA_n(y,ff))) && CAA_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_CAA(ff) -= comp_multinom(CAA_hist, CAApred, CN, CAA_n, y, max_age, ff);
          } else {
            nll_CAA(ff) -= comp_lognorm(CAA_hist, CAApred, CN, CAA_n, y, max_age, ff);
          }
        }

        if(!R_IsNA(asDouble(CAL_n(y,ff))) && CAL_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_CAL(ff) -= comp_multinom(CAL_hist, CALpred, CN, CAL_n, y, nlbin, ff);
          } else {
            nll_CAL(ff) = comp_lognorm(CAL_hist, CALpred, CN, CAL_n, y, nlbin, ff);
          }
        }

        if(condition == "catch" || nll_C) nll_Catch(ff) -= dnorm(log(C_hist(y,ff)), log(Cpred(y,ff)), Type(0.01), true);
        if(!R_IsNA(asDouble(mlen(y,ff))) && mlen(y,ff) > 0) {
          nll_ML(ff) -= dnorm(mlen(y,ff), mlen_pred(y,ff), sigma_mlen(ff), true);
        }
      }
    }

    nll_Catch(ff) *= LWT_C(ff,0);
    nll_CAA(ff) *= LWT_C(ff,1);
    nll_CAL(ff) *= LWT_C(ff,2);
    nll_ML(ff) *= LWT_C(ff,3);

    if(condition == "catch" && C_eq(ff) > 0) nll_Ceq(ff) = -1 * LWT_C(ff,4) * dnorm(log(C_eq(ff)), log(C_eq_pred(ff)), Type(0.01), true);

  }

  for(int y=0;y<n_y;y++) {
    if(est_rec_dev(y)) nll_log_rec_dev -= dnorm(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<max_age-1;a++) {
    if(est_early_rec_dev(a)) nll_log_rec_dev -= dnorm(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_Catch.sum() + nll_Index.sum();
  nll += nll_s_CAA.sum() + nll_s_CAL.sum();
  nll += nll_CAA.sum() + nll_CAL.sum() + nll_ML.sum();
  nll += nll_log_rec_dev + nll_Ceq.sum();
  nll += penalty + prior;

  if(CppAD::Variable(log_R0)) ADREPORT(R0);
  ADREPORT(h);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(condition == "effort") ADREPORT(q_effort);
  if(nll_Index.sum() != 0) ADREPORT(q);

  REPORT(M);
  REPORT(length_bin);
  REPORT(Linf);

  REPORT(log_R0);
  REPORT(transformed_h);
  REPORT(vul_par);
  REPORT(LFS);
  REPORT(L5);
  REPORT(Vmaxlen);
  REPORT(log_q_effort);
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
  REPORT(Z);

  REPORT(NPR_unfished);
  REPORT(EPR0);
  REPORT(B0);
  REPORT(E0);
  REPORT(N0);

  REPORT(Arec);
  REPORT(Brec);
  REPORT(E0_SR);
  REPORT(EPR0_SR);

  REPORT(ALK);
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
  REPORT(EPR_eq);
  REPORT(R_eq);

  REPORT(C_eq_pred);
  REPORT(q);

  REPORT(nll_Catch);
  REPORT(nll_Index);
  REPORT(nll_s_CAA);
  REPORT(nll_s_CAL);
  REPORT(nll_CAA);
  REPORT(nll_CAL);
  REPORT(nll_ML);
  REPORT(nll_Ceq);
  REPORT(nll_log_rec_dev);

  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  REPORT(s_CAApred);
  REPORT(s_CALpred);
  REPORT(s_vul);
  REPORT(s_L5);
  REPORT(s_LFS);
  REPORT(s_Vmaxlen);

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif



