
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace density;

  DATA_SCALAR(S0);
  DATA_SCALAR(Alpha);
  DATA_SCALAR(Rho);
  DATA_INTEGER(ny);
  DATA_INTEGER(k);
  DATA_SCALAR(wk);
  DATA_VECTOR(E_hist);
  DATA_VECTOR(C_hist);
  DATA_STRING(SR_type);
  DATA_VECTOR_INDICATOR(keep, C_hist);

  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(log_q);
  PARAMETER(U_equilibrium);
  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_rec_dev);

  Type UMSY = ilogit(logit_UMSY);
  Type MSY = exp(log_MSY);
  Type q = exp(log_q);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);

  //--DECLARING DERIVED VALUES
  Type BMSY = MSY/UMSY;
  Type SS = S0 * (1 - UMSY);
  Type Spr = (SS * Alpha/(1 - SS) + wk)/(1 - Rho * SS);
  Type DsprDu = (Alpha + Spr * (1 + Rho - 2 * Rho * SS))/((1 - Rho * SS) * (1 - SS));
  DsprDu += Alpha * SS/((1 - Rho * SS) * pow(1 - SS, 2));
  DsprDu -= Spr/(1 - SS);
  DsprDu *= -S0;
  Type Spr0 = (S0 * Alpha/(1 - S0) + wk)/(1 - Rho * S0);

  Type Arec;
  Type Brec;
  Type h;
  Type R0;

  if(SR_type == "BH") {
    Arec = 1/(pow(1 - UMSY, 2) * (Spr + UMSY * DsprDu));
    Brec = UMSY * (Arec * Spr - 1/(1 - UMSY))/MSY;
    h = Arec * Spr0 / (4 + Arec * Spr0);
    R0 = (Arec * Spr0 - 1)/(Brec * Spr0);
  } else {
    Type log_aphi = -UMSY * (1 - UMSY) * DsprDu;
    log_aphi /=	Spr;
    log_aphi += UMSY;

    Arec = exp(log_aphi);
    Arec /= Spr * (1 - UMSY);
    Brec = UMSY * Spr * log_aphi;
    Brec /= MSY * (1 - UMSY);

    h = exp(0.8 * Brec * Spr0);
    h *= 0.2;
    R0 = log(Arec * Spr0)/(Brec * Spr0);
  }

  Type B0 = R0 * Spr0;
  Type N0 = R0/(1 - S0);

  //--DECLARING STORAGE VECTORS
  int ny_p = ny + 1;
  int ny_k = ny + k;
  vector<Type> B(ny_p);
  vector<Type> N(ny_p);
  vector<Type> R(ny_k);
  vector<Type> Rec_dev(ny - k);

  vector<Type> Surv(ny);
  vector<Type> Cpred(ny);
  vector<Type> Sp(ny);
  vector<Type> U(ny);

  //--INITIALIZE
  Type Seq = S0 * (1 - U_equilibrium);
  Type SprEq = (Seq * Alpha/(1 - Seq) + wk)/(1 - Rho * Seq);
  Type Req;
  if(SR_type == "BH") {
    Req = (Arec * SprEq * (1 - U_equilibrium) - 1)/(Brec * SprEq * (1 - U_equilibrium));
  } else {
    Req = log(Arec * SprEq * (1 - U_equilibrium))/(Brec * SprEq * (1 - U_equilibrium));
  }

  B(0) = Req * SprEq;
  N(0) = Req/(1 - Seq);
  for(int tt=0;tt<k;tt++) R(tt) = Req;

  Type penalty = 0;

  for(int tt=0; tt<ny; tt++){
    U(tt) = CppAD::CondExpLt(exp(-q * E_hist(tt)), Type(0.025),
      1 - posfun(exp(-q * E_hist(tt)), Type(0.025), penalty), 1 - exp(-q * E_hist(tt)));
    Surv(tt) = S0 * (1 - U(tt));
    Cpred(tt) = U(tt) * B(tt);
    Sp(tt) = B(tt) - Cpred(tt);

    if(SR_type == "BH") {
      R(tt + k) = Arec * Sp(tt)/(1 + Brec * Sp(tt));
    } else {
      R(tt + k) = Arec * Sp(tt) * exp(-Brec * Sp(tt));
    }
    if(tt + k < ny) {
      Rec_dev(tt) = exp(log_rec_dev(tt) - 0.5 * pow(tau, 2));
      R(tt + k) *= Rec_dev(tt);
    }

    B(tt+1) = Surv(tt) * (Alpha * N(tt) + Rho * B(tt)) + wk * R(tt+1);
    N(tt+1) = Surv(tt) * N(tt) + R(tt+1);
  }

  //--ARGUMENTS FOR NLL
  // The following conditions must be met for positive values of:
  // Arec: UMSY * DsprDu + Spr > 0 (BH)
  //       Spr * (1 - UMSY) > 0 (Ricker)
  // Brec: Arec * Spr * (1 - UMSY) - 1 > 0 (both BH and Ricker)
  // Thus, create a likelihood penalty of 100 when either condition is not met
  penalty += CppAD::CondExpGt(Arec * Spr * (1 - UMSY) - 1, Type(0), Type(0), Type(UMSY * 1e3));
  if(SR_type == "BH") penalty += CppAD::CondExpGt(Spr + UMSY * DsprDu, Type(0), Type(0), Type(UMSY * 1e3));
  if(SR_type == "Ricker") penalty += CppAD::CondExpGt(Spr * (1 - UMSY), Type(0), Type(0), Type(UMSY * 1e3));

  // Objective function
  //creates storage for jnll and sets value to 0
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();

  for(int tt=0; tt<ny; tt++){
    jnll_comp(0) -= keep(tt) * dnorm(log(C_hist(tt)), log(Cpred(tt)), sigma, true);
    if(tt+k<ny) jnll_comp(1) -= dnorm(log_rec_dev(tt), Type(0), tau, true);
  }

  //Summing individual jnll and penalties
  Type jnll = jnll_comp.sum() + penalty;

  //-------REPORTING-------//
  Type U_UMSY_final = U(U.size()-1)/UMSY;
  Type B_BMSY_final = B(B.size()-1)/BMSY;
  Type B_B0_final = B(B.size()-1)/B0;

  ADREPORT(UMSY);
  ADREPORT(MSY);
  ADREPORT(q);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(U_UMSY_final);
  ADREPORT(B_BMSY_final);
  ADREPORT(B_B0_final);
  REPORT(UMSY);
  REPORT(MSY);
  REPORT(sigma);
  REPORT(tau);
  REPORT(jnll_comp);
  REPORT(jnll);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(Spr);
  REPORT(Spr0);
  REPORT(DsprDu);
  REPORT(Cpred);
  REPORT(B);
  REPORT(N);
  REPORT(R);
  REPORT(log_rec_dev);
  REPORT(Rec_dev);
  REPORT(U);
  REPORT(h);
  REPORT(BMSY);
  REPORT(R0);
  REPORT(N0);
  REPORT(B0);
  REPORT(penalty);

  return jnll;
//}
