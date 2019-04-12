
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

  PARAMETER(log_R0);
  PARAMETER(transformed_h);
  PARAMETER(log_q);
  PARAMETER(U_equilibrium);

  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else h = exp(transformed_h);
  h += 0.2;
  Type R0 = exp(log_R0);
  Type q = exp(log_q);

  //--DECLARING DERIVED VALUES
  Type Spr0 = (S0 * Alpha/(1 - S0) + wk)/(1 - Rho * S0);
  Type B0 = R0 * Spr0;
  Type N0 = R0/(1 - S0);

  Type Arec;
  Type Brec;

  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Arec /= Spr0;
    Brec = 5*h - 1;
    Brec /= (1-h) * B0;
  } else {
    Arec = pow(5*h, 1.25);
    Arec /= Spr0;
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= B0;
  }

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
      R(tt + k) = Arec * B(tt)/(1 + Brec * B(tt));
    } else {
      R(tt + k) = Arec * B(tt) * exp(-Brec * B(tt));
    }
    B(tt + 1) = Surv(tt) * (Alpha * N(tt) + Rho * B(tt)) + wk * R(tt + 1);
    N(tt + 1) = Surv(tt) * N(tt) + R(tt + 1);
  }

  //--ARGUMENTS FOR NLL
  Type sigma = calc_sigma(C_hist, Cpred);

  // Objective function
  //creates storage for jnll and sets value to 0
  Type nll = 0.;
  for(int tt=0; tt<ny; tt++) {
    if(C_hist(tt) > 0) nll -= dnorm(log(C_hist(tt)), log(Cpred(tt)), sigma, true);
  }

  //Summing individual jnll and penalties
  nll += penalty;

  //-------REPORTING-------//
  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(q);
  ADREPORT(sigma);
  REPORT(sigma);
  REPORT(nll);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(Spr0);
  REPORT(Cpred);
  REPORT(B);
  REPORT(N);
  REPORT(R);
  REPORT(U);
  REPORT(h);
  REPORT(R0);
  REPORT(N0);
  REPORT(B0);
  REPORT(penalty);

  return nll;
//}
