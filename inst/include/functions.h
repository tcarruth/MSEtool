// Calculates analytical solution of a lognormal variable
template<class Type>
Type calc_sigma(vector<Type> I_y, vector<Type> Ipred_y) {
  Type sum_square = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.size();y++) {
    if(I_y(y)>0) {
      sum_square += pow(log(I_y(y)/Ipred_y(y)), 2);
      n_y += 1.;
    }
  }
  Type sigma = pow(sum_square/n_y, 0.5);
  return sigma;
}

// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(vector<Type> I_y, vector<Type> B_y) {
  Type num = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.size();y++) {
    if(I_y(y)>0) {
      num += log(I_y(y)/B_y(y));
      n_y += 1.;
    }
  }
  Type q = exp(num/n_y);
  return q;
}

