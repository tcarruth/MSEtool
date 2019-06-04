
namespace SRA_scope {

using namespace SCA;

template<class Type>
matrix<Type> generate_ALK(vector<Type> length_bin, matrix<Type> len_age, Type CV_LAA, int max_age, int nlbin, Type bin_width, int y) {
  matrix<Type> ALK(max_age, length_bin.size());
  for(int a=0;a<max_age;a++) {
    for(int j=0;j<nlbin;j++) {
      if(j==nlbin-1) {
        ALK(a,j) = 1 - pnorm(length_bin(j) - 0.5 * bin_width, len_age(y,a), CV_LAA * len_age(y,a));
      } else {
        ALK(a,j) = pnorm(length_bin(j) + 0.5 * bin_width, len_age(y,a), CV_LAA * len_age(y,a));
        if(j>0) ALK(a,j) -= pnorm(length_bin(j) - 0.5 * bin_width, len_age(y,a), CV_LAA * len_age(y,a));
      }
    }
  }
  return ALK;
}

template<class Type>
matrix<Type> calc_NPR0(int nlbin, matrix<Type> M, int max_age, matrix<Type> ALK, int y) {
  vector<Type> NPR(max_age);
  NPR.setZero();
  matrix<Type> NPR_full(max_age, nlbin);

  NPR(0) = 1;
  for(int len=0;len<nlbin;len++) NPR_full(0,len) = NPR(0) * ALK(0,len);
  for(int a=1;a<max_age;a++) {
    for(int len=0;len<nlbin;len++) NPR(a) += NPR_full(a-1,len) * exp(-M(y,a-1));
    for(int len=0;len<nlbin;len++) NPR_full(a,len) = NPR(a) * ALK(a,len);
  }
  for(int len=0;len<nlbin;len++) {
    NPR_full(max_age-1,len) /= 1 - exp(-M(y,max_age-1));
    NPR(max_age-1) += NPR_full(max_age-1,len);
  }
  return NPR_full;
}


template<class Type>
matrix<Type> calc_NPR(vector<Type> F, vector<vector<Type> > vul, int nfleet, int nlbin, matrix<Type> M, int max_age, matrix<Type> ALK, int y) {
  vector<Type> NPR(max_age);
  NPR.setZero();
  matrix<Type> NPR_full(max_age, nlbin);

  NPR(0) = 1;
  for(int len=0;len<nlbin;len++) NPR_full(0,len) = NPR(0) * ALK(0,len);
  for(int a=1;a<max_age;a++) {
    for(int len=0;len<nlbin;len++) {
      Type Z_total = M(y,a-1);
      for(int ff=0;ff<nfleet;ff++) Z_total += vul(ff)(len) * F(ff);
      NPR(a) += NPR_full(a-1,len) * exp(-Z_total);
    }
    for(int len=0;len<nlbin;len++) NPR_full(a,len) = NPR(a) * ALK(a,len);
  }

  for(int len=0;len<nlbin;len++) {
    Type Z_total = M(y,max_age-1);
    for(int ff=0;ff<nfleet;ff++) Z_total += vul(ff)(len) * F(ff);
    NPR_full(max_age-1,len) /= 1 - exp(-Z_total);
    NPR(max_age-1) += NPR_full(max_age-1,len);
  }
  return NPR_full;
}

template<class Type>
Type get_N_at_age(matrix<Type> NPR, int nlbin, int a) {
  Type N = 0;
  for(int len=0;len<nlbin;len++) N += NPR(a,len);
  return N;
}

template<class Type>
Type sum_EPR(matrix<Type> NPR, vector<Type> wt_at_len, matrix<Type> mat, int max_age, int nlbin, int y) {
  Type EPR = 0.;
  for(int a=0;a<max_age;a++) {
    for(int len=0;len<nlbin;len++) EPR += NPR(a,len) * wt_at_len(len) * mat(y,a);
  }
  return EPR;
}

template<class Type>
Type sum_BPR(matrix<Type> NPR, vector<Type> wt_at_len) {
  Type BPR = (NPR * wt_at_len).sum();
  return BPR;
}

template<class Type>
Type sum_VBPR(matrix<Type> NPR, vector<Type> wt_at_len, vector<Type> vul, int max_age, int nlbin) {
  Type VBPR = 0;
  for(int a=0;a<max_age;a++) {
    for(int len=0;len<nlbin;len++) VBPR += NPR(a,len) * wt_at_len(len) * vul(len);
  }
  return VBPR;
}

template<class Type>
vector<Type> calc_C_eq(vector<Type> F, array<Type> N, vector<vector<Type> > vul, matrix<Type> M, vector<Type> wt_at_len, int nlbin,
                       int nfleet, int max_age, int y) {
  vector<Type> C_eq(nfleet);
  C_eq.setZero();
  for(int a=0;a<max_age;a++) {
    for(int len=0;len<nlbin;len++) {
      Type Z_total = 0;
      for(int ff=0;ff<nfleet;ff++) Z_total += vul(ff)(len) * F(ff) + M(y,a);
      for(int ff=0;ff<nfleet;ff++) C_eq(ff) += vul(ff)(len) * F(ff) * N(y,a,len) * (1 - exp(-Z_total))/Z_total * wt_at_len(len);
    }
  }
  return C_eq;
}





template<class Type>
vector<Type> calc_dome_vul(vector<Type> vul_par, int nlbin, vector<Type> length_bin, Type &prior) {
  vector<Type> vul(nlbin);

  Type max_lenbin = length_bin(length_bin.size()-1);

  Type len_full = invlogit(vul_par(0)) * 0.75 * max_lenbin;
  Type len_50 = len_full - exp(vul_par(1));
  Type len_full2 = invlogit(vul_par(2));
  len_full2 *= max_lenbin - len_full;
  len_full2 += len_full;
  Type vul_max = invlogit(vul_par(3));

  //prior -= dnorm(vul_par(0), Type(0), Type(3), true);
  prior -= dnorm(vul_par(1), Type(0), Type(3), true);
  prior -= dnorm(vul_par(2), Type(0), Type(3), true);
  //prior -= dnorm(vul_par(3), Type(0), Type(3), true);

  Type var_asc = pow(len_50 - len_full, 2);
  var_asc /= log(Type(4));

  Type var_des = pow(max_lenbin - len_full2, 2);
  var_des /= -2 * log(vul_max);

  Type sd_asc = pow(var_asc, 0.5);
  Type sd_des = pow(var_des, 0.5);

  for(int len=0;len<nlbin;len++) {
    Type vul_asc = dnorm_vul(length_bin(len), len_full, sd_asc);
    Type vul_des = dnorm_vul(length_bin(len), len_full2, sd_des);

    vul(len) = CppAD::CondExpLe(length_bin(len), len_full, vul_asc, CppAD::CondExpLe(length_bin(len), len_full2, Type(1), vul_des));
  }
  Type interim_vmax = max(vul);
  for(int len=0;len<nlbin;len++) vul(len) /= interim_vmax;
  return vul;
}



template<class Type>
vector<Type> calc_logistic_vul(vector<Type> vul_par, int nlbin, vector<Type> length_bin, Type &prior) {
  vector<Type> vul(nlbin);
  Type max_lenbin = length_bin(length_bin.size()-1);
  Type len_95 = invlogit(vul_par(0)) * 0.75 * max_lenbin;
  Type len_50 = len_95 - exp(vul_par(1));

  prior -= dnorm(vul_par(1), Type(0), Type(3), true);

  for(int len=0;len<nlbin;len++) {
    vul(len) = pow(1 + exp(-log(Type(19.0)) * (length_bin(len) - len_50)/(len_95 - len_50)), -1);
  }
  Type interim_vmax = max(vul);
  for(int len=0;len<nlbin;len++) vul(len) /= interim_vmax;
  return vul;
}

template<class Type>
vector<vector<Type> > calc_vul(matrix<Type> vul_par, vector<int> vul_type, int nlbin, vector<Type> length_bin, Type &prior) {
  vector<vector<Type> > vul(vul_type.size());

  for(int ff=0;ff<vul_type.size();ff++) {
    vector<Type> vul_par_ff(4);
    vul_par_ff(0) = vul_par(0,ff);
    vul_par_ff(1) = vul_par(1,ff);
    vul_par_ff(2) = vul_par(2,ff);
    vul_par_ff(3) = vul_par(3,ff);
    if(vul_type(ff)) {
      vul(ff) = calc_logistic_vul(vul_par_ff, nlbin, length_bin, prior);
    } else {
      vul(ff) = calc_dome_vul(vul_par_ff, nlbin, length_bin, prior);
    }
  }
  return vul;
}


// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(matrix<Type> I_y, vector<Type> B_y, int sur, matrix<Type> &Ipred) {
  Type num = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.rows();y++) {
    if(!R_IsNA(asDouble(I_y(y,sur))) && I_y(y,sur)>0) {
      num += log(I_y(y,sur)/B_y(y));
      n_y += 1.;
    }
  }
  Type q = exp(num/n_y);
  for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q * B_y(y);
  return q;
}


// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(matrix<Type> I_y, matrix<Type> B_y, int sur, int ff, matrix<Type> &Ipred) {
  Type num = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.rows();y++) {
    if(!R_IsNA(asDouble(I_y(y,sur))) && I_y(y,sur)>0) {
      num += log(I_y(y,sur)/B_y(y,ff));
      n_y += 1.;
    }
  }
  Type q = exp(num/n_y);
  for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q * B_y(y,ff);
  return q;
}


}
