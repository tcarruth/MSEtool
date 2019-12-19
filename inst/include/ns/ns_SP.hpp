
namespace ns_SP {

// Sub-annual time steps for SP and iteratively find the F that predicts the observed catch
template<class Type>
Type SP_F(Type U_start, Type C_hist, Type MSY, Type K, Type n, Type nterm, Type dt, int nstep, int nitF,
          vector<Type> &Cpred, vector<Type> &B, int y, Type &penalty) {

  Type F = -log(1 - U_start);
  if(dt < 1) {
    for(int i=0;i<nitF;i++) {
      Type Catch = 0;
      Type B_next = B(y);

      for(int seas=0;seas<nstep;seas++) {
        Type SP = CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B_next / K * log(B_next/K),
                                   nterm/(n-1) * MSY * (B_next/K - pow(B_next/K, n))) - F * B_next;
        SP *= dt;
        Catch += F * B_next * dt;
        B_next += SP;
      }

      if(i==nitF-1) {
        B(y+1) = B_next;
        Cpred(y) = Catch;
      } else {
        F *= C_hist/Catch;
      }
    }
  } else {
    Cpred(y) = C_hist;
    F = CppAD::CondExpLt(1 - C_hist/B(y), Type(0.025), 1 - posfun(1 - C_hist/B(y), Type(0.025), penalty), C_hist/B(y));
    B(y+1) = B(y) + CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B(y) / K * log(B(y)/K),
      nterm/(n-1) * MSY * (B(y)/K - pow(B(y)/K, n))) - F * B(y);
  }

  return F;
}

template<class Type>
vector<Type> calc_q(matrix<Type> I_y, vector<Type> B_y, matrix<Type> &Ipred, int nsurvey) {
  vector<Type> num(nsurvey);
  vector<Type> n_y(nsurvey);
  num.setZero();
  n_y.setZero();

  vector<Type> q(nsurvey);

  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<I_y.rows();y++) {
      if(!R_IsNA(asDouble(I_y(y,sur))) && I_y(y,sur)>0) {
        num(sur) += log(I_y(y,sur)/B_y(y));
        n_y(sur) += 1.;
      }
    }
    q(sur) = exp(num(sur)/n_y(sur));
    for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q(sur) * B_y(y);
  }
  return q;
}


}
