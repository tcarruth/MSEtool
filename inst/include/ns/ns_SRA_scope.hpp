
namespace ns_SRA_scope {

using namespace ns_SCA;

template <class Type>
Type log2(Type x) {
  return log(x)/log(Type(2));
}

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
vector<Type> calc_NPR0(matrix<Type> M, int max_age, int y) {
  vector<Type> NPR(max_age);
  NPR(0) = 1;
  for(int a=1;a<max_age;a++) NPR(a) = NPR(a-1) * exp(-M(y,a-1));
  NPR(max_age-1) /= 1 - exp(-M(y,max_age-1));
  return NPR;
}


template<class Type>
vector<Type> calc_NPR(vector<Type> F, array<Type> vul, int nfleet, matrix<Type> M, int max_age, int y) {
  vector<Type> NPR(max_age);
  NPR(0) = 1;
  for(int a=1;a<max_age;a++) {
    Type Z = M(y,a);
    for(int ff=0;ff<nfleet;ff++) Z += vul(y,a,ff) * F(ff);
    NPR(a) = NPR(a-1) * exp(-Z);
    if(a == max_age-1) NPR(a) /= 1 - exp(-Z);
  }
  return NPR;
}

template<class Type>
Type sum_EPR(vector<Type> NPR, matrix<Type> wt, matrix<Type> mat, int max_age, int y) {
  Type EPR = 0.;
  for(int a=0;a<max_age;a++) EPR += NPR(a) * wt(y,a) * mat(y,a);
  return EPR;
}

template<class Type>
Type sum_BPR(vector<Type> NPR, matrix<Type> wt, int max_age, int y) {
  Type BPR = 0;
  for(int a=0;a<max_age;a++) BPR += NPR(a) * wt(y,a);
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
array<Type> calc_vul(matrix<Type> vul_par, vector<int> vul_type, matrix<Type> Len_age, vector<Type> &LFS, vector<Type> &L5,
                     vector<Type> &Vmaxlen, Type Linf) {
  array<Type> vul(Len_age.rows(), Len_age.cols(), vul_type.size());

  for(int ff=0;ff<vul_type.size();ff++) {
    LFS(ff) = invlogit(vul_par(0,ff)) * 0.95 * Linf;
    L5(ff) = LFS(ff) - exp(vul_par(1,ff));
    Type sls = (LFS(ff) - L5(ff))/pow(-log2(0.05), 0.5);

    if(vul_type(ff) < 0) { // Logistic
      Vmaxlen(ff) = 1;

      for(int y=0;y<Len_age.rows();y++) {
        for(int a=0;a<Len_age.cols();a++) {
          Type lo = pow(2, -((Len_age(y,a) - LFS(ff))/sls * (Len_age(y,a) - LFS(ff))/sls));
          vul(y,a,ff) = CppAD::CondExpLt(Len_age(y,a), LFS(ff), lo, Type(1));
        }
      }
    } else if(vul_type(ff) == 0) { // Dome
      Vmaxlen(ff) = invlogit(vul_par(2,ff));
      Type srs = (Linf - LFS(ff))/pow(-log2(Vmaxlen(ff)), 0.5);

      for(int y=0;y<Len_age.rows();y++) {
        for(int a=0;a<Len_age.cols();a++) {
          Type lo = pow(2, -((Len_age(y,a) - LFS(ff))/sls * (Len_age(y,a) - LFS(ff))/sls));
          Type hi = pow(2, -((Len_age(y,a) - LFS(ff))/srs * (Len_age(y,a) - LFS(ff))/srs));
          vul(y,a,ff) = CppAD::CondExpLt(Len_age(y,a), LFS(ff), lo, hi);
        }
      }
    } else { // Age-specific index
      for(int y=0;y<Len_age.rows();y++) {
        for(int a=0;a<Len_age.cols();a++) {
          if(a == vul_type(ff) - 1) {
            vul(y,a,ff) = 1;
          } else {
            vul(y,a,ff) = 0;
          }
        }
      }
    }
  }

  return vul;
}


// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(matrix<Type> I_y, vector<Type> B_y, int sur, matrix<Type> &Ipred, vector<int> abs_I) {

  Type q;
  if(abs_I(sur)) {
    q = 1;
  } else {
    Type num = 0.;
    Type n_y = 0.;

    for(int y=0;y<I_y.rows();y++) {
      if(!R_IsNA(asDouble(I_y(y,sur))) && I_y(y,sur)>0) {
        num += log(I_y(y,sur)/B_y(y));
        n_y += 1.;
      }
    }
    q = exp(num/n_y);
  }
  for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q * B_y(y);
  return q;
}


// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(matrix<Type> I_y, matrix<Type> B_y, int sur, int ff, matrix<Type> &Ipred, vector<int> abs_I) {
  Type q;
  if(abs_I(sur)) {
    q = 1;
  } else {
    Type num = 0.;
    Type n_y = 0.;

    for(int y=0;y<I_y.rows();y++) {
      if(!R_IsNA(asDouble(I_y(y,sur))) && I_y(y,sur)>0) {
        num += log(I_y(y,sur)/B_y(y,ff));
        n_y += 1.;
      }
    }
    q = exp(num/n_y);
  }
  for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q * B_y(y,ff);
  return q;
}


}
