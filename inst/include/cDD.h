
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace density;

  DATA_SCALAR(M);
  DATA_SCALAR(Winf);
  DATA_SCALAR(Kappa);
  DATA_INTEGER(ny);
  DATA_INTEGER(k);
  DATA_SCALAR(wk);
  DATA_VECTOR(C_hist);
  DATA_VECTOR(I_hist);
  DATA_STRING(SR_type);
  DATA_INTEGER(nitF);

  PARAMETER(log_R0);
  PARAMETER(transformed_h);
  PARAMETER(F_equilibrium);

  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else h = exp(transformed_h);
  h += 0.2;
  Type R0 = exp(log_R0);
  int SR_type2 = SR_type == "BH";

  //--DECLARING DERIVED VALUES
  Type BPR0 = cDD_BPR(Type(0), M, wk, Kappa, Winf);
  Type B0 = BPR0 * R0;
  Type N0 = R0/M;

  Type Arec;
  Type Brec;

  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Arec /= BPR0;
    Brec = 5*h - 1;
    Brec /= (1-h) * B0;
  } else {
    Arec = pow(5*h, 1.25);
    Arec /= BPR0;
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= B0;
  }

  //--DECLARING STORAGE VECTORS
  vector<Type> B(ny+1);
  vector<Type> N(ny+1);
  vector<Type> R(ny+k);

  vector<Type> F(ny);
  vector<Type> Z(ny);
  vector<Type> Cpred(ny);
  vector<Type> Ipred(ny);

  vector<Type> BPRinf(ny);
  vector<Type> Rinf(ny);
  vector<Type> Binf(ny);
  vector<Type> Ninf(ny);

  //--INITIALIZE
  Type BPReq = cDD_BPR(F_equilibrium, M, wk, Kappa, Winf);
  Type Req = cDD_R(BPReq, Arec, Brec, SR_type2);

  B(0) = Req * BPReq;
  N(0) = Req/(F_equilibrium + M);
  for(int tt=0;tt<k;tt++) R(tt) = Req;

  Type penalty = 0;

  for(int tt=0; tt<ny; tt++) {
	F(tt) = cDD_F(C_hist(tt)/B(tt), C_hist(tt), M, Winf, Kappa, wk, Arec, Brec, SR_type2, N, B, Cpred, BPRinf, Binf, R, Ninf, nitF, tt);
	Z(tt) = F(tt) + M;

	N(tt+1) = Ninf(tt) + (N(tt) - Ninf(tt)) * exp(-Z(tt));
	B(tt+1) = Binf(tt);
	B(tt+1) += Kappa * Winf * (N(tt) - Ninf(tt)) / (Z(tt) + Kappa);
	B(tt+1) += (B(tt) - Binf(tt) - Kappa * Winf * (N(tt) - Ninf(tt))/(Z(tt) + Kappa)) * exp(-Z(tt) - Kappa);

	if(SR_type == "BH") {
      R(tt+k) = Arec * B(tt)/(1 + Brec * B(tt));
    } else {
      R(tt+k) = Arec * B(tt) * exp(-Brec * B(tt));
    }
  }

  //--ARGUMENTS FOR NLL
  Type q = calc_q(I_hist, B);
  for(int y=0;y<ny;y++) Ipred(y) = q * B(y);
  Type sigma = calc_sigma(I_hist, Ipred);

  // Objective function
  //creates storage for jnll and sets value to 0
  Type nll = 0.;
  for(int y=0;y<ny;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
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
  REPORT(BPR0);
  REPORT(Cpred);
  REPORT(Ipred);
  REPORT(B);
  REPORT(N);
  REPORT(R);
  REPORT(F);
  REPORT(Z);
  REPORT(BPRinf);
  REPORT(Binf);
  REPORT(Ninf);
  REPORT(Rinf);
  REPORT(h);
  REPORT(R0);
  REPORT(N0);
  REPORT(B0);
  REPORT(penalty);

  return nll;
//}
