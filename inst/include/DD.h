
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

  PARAMETER(logit_UMSY);
  PARAMETER(log_MSY);
  PARAMETER(log_q);

  Type UMSY = 1/(1 + exp(-logit_UMSY));
  Type MSY = exp(log_MSY);
  Type q = exp(log_q);

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

  vector<Type> Surv(ny);
  vector<Type> Cpred(ny);
  vector<Type> Sp(ny);
  vector<Type> U(ny);

  //--INITIALIZE
  B(0) = B0;
  N(0) = N0;
  for(int tt=0;tt<k;tt++) R(tt) = R0;

  for(int tt=0; tt<ny; tt++){
    U(tt) = 1 - exp(-q * E_hist(tt));
    Surv(tt) = S0 * (1 - U(tt));
    Cpred(tt) = CppAD::CondExpGt(U(tt) * B(tt), Type(1e-15), U(tt) * B(tt), Type(1e-15));
    Sp(tt) = B(tt) - Cpred(tt);

    if(SR_type == "BH") {
      R(tt + k) = Arec * Sp(tt)/(1 + Brec * Sp(tt));
    } else {
      R(tt + k) = Arec * Sp(tt) * exp(-Brec * Sp(tt));
    }
    B(tt + 1) = Surv(tt) * (Alpha * N(tt) + Rho * B(tt)) + wk * R(tt + 1);
    N(tt + 1) = Surv(tt) * N(tt) + R(tt + 1);
  }

  //--ARGUMENTS FOR NLL
  Type sigma = calc_sigma(C_hist, Cpred);

  // The following conditions must be met for positive values of:
  // Arec: UMSY * DsprDu + Spr > 0 (BH)
  //       Spr * (1 - UMSY) > 0 (Ricker)
  // Brec: Arec * Spr * (1 - UMSY) - 1 > 0 (both BH and Ricker)
  // Thus, create a likelihood penalty of 100 when either condition is not met
  Type penalty = CppAD::CondExpGt(Arec * Spr * (1 - UMSY) - 1, Type(0), Type(0), Type(UMSY * 1e3));
  if(SR_type == "BH") penalty += CppAD::CondExpGt(Spr + UMSY * DsprDu, Type(0), Type(0), Type(UMSY * 1e3));
  if(SR_type == "Ricker") penalty += CppAD::CondExpGt(Spr * (1 - UMSY), Type(0), Type(0), Type(UMSY * 1e3));

  // Objective function
  //creates storage for jnll and sets value to 0
  Type nll = 0.;
  for(int tt=0; tt<ny; tt++) nll -= dnorm(log(C_hist(tt)), log(Cpred(tt)), sigma, true);

  //Summing individual jnll and penalties
  nll += penalty;

  //-------REPORTING-------//
  Type U_UMSY_final = U(U.size()-1)/UMSY;
  Type B_BMSY_final = B(B.size()-1)/BMSY;
  Type B_B0_final = B(B.size()-1)/B0;

  ADREPORT( UMSY );
  ADREPORT( MSY );
  ADREPORT( q );
  ADREPORT( sigma );
  ADREPORT( U_UMSY_final );
  ADREPORT( B_BMSY_final );
  ADREPORT( B_B0_final );
  REPORT( UMSY );
  REPORT( MSY );
  REPORT( q );
  REPORT( sigma );
  REPORT( nll );
  REPORT( Arec );
  REPORT( Brec );
  REPORT( Spr );
  REPORT( Spr0 );
  REPORT( DsprDu );
  REPORT( Cpred );
  REPORT( B );
  REPORT( N );
  REPORT( R );
  REPORT( U );
  REPORT( h );
  REPORT( BMSY );
  REPORT( R0 );
  REPORT( N0 );
  REPORT( B0 );
  REPORT( penalty );

  return nll;
//}
