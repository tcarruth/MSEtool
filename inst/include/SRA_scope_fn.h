

template<class Type>
Type calc_q(matrix<Type> I_y, vector<Type> B_y, int ff) {
  
  vector<Type> I_hist(B_y.size()-1);
  for(int y=0;y<I_hist.size();y++) I_hist(y) = I_y(y,ff);
  
  return calc_q(I_hist, B_y);
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
matrix<Type> calc_NPR(Type U, vector<Type> vul, int nlbin, matrix<Type> M, int max_age, matrix<Type> ALK, int y) {
  vector<Type> NPR(max_age);
  NPR.setZero();
  matrix<Type> NPR_full(max_age, nlbin);
  
  NPR(0) = 1;
  for(int len=0;len<nlbin;len++) NPR_full(0,len) = NPR(0) * ALK(0,len);
  for(int a=1;a<max_age;a++) {
	for(int len=0;len<nlbin;len++) NPR(a) += NPR_full(a-1,len) * exp(-M(y,a-1)) * (1 - vul(len) * U);
	for(int len=0;len<nlbin;len++) NPR_full(a,len) = NPR(a) * ALK(a,len);
  }
  for(int len=0;len<nlbin;len++) {
	NPR_full(max_age-1,len) /= 1 - exp(-M(y,max_age-1)) * (1 - vul(len) * U);
	NPR(max_age-1) += NPR_full(max_age-1,len);
  }
  return NPR_full;  
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
  //Type answer = 0.;
  //for(int a=0;a<NPR.rows();a++) answer += NPR(a) * weight(y,a);
  //return answer;
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
  
  //vector<Type> answer(nlbin);
  //answer = NPR * wt_at_len;
  //return (answer * vul).sum();
}

template<class Type>
Type calc_C_eq(Type U, array<Type> N, vector<Type> vul, matrix<Type> M, vector<Type> wt_at_len, int nlbin, int max_age, int y) {
  Type C_eq = 0;
  for(int a=0;a<max_age;a++) {
    for(int len=0;len<nlbin;len++) C_eq += U * N(y,a,len) * exp(0.5 * M(y,a)) * vul(len) * wt_at_len(len);
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
