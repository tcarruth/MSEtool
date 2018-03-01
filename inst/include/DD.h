
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace density;

  DATA_SCALAR( So_DD );
  DATA_SCALAR( Alpha_DD );
  DATA_SCALAR( Rho_DD );
  DATA_INTEGER( ny_DD );
  DATA_INTEGER( k_DD );
  DATA_SCALAR( wa_DD );
  DATA_VECTOR( E_hist );
  DATA_VECTOR( C_hist );
  //DATA_VECTOR( UMSYprior );

  PARAMETER( logit_UMSY_DD );
  PARAMETER( log_MSY_DD );
  PARAMETER( log_q_DD );

  Type UMSY_DD = 1/(1 + exp(-logit_UMSY_DD));
  Type MSY_DD = exp(log_MSY_DD);
  Type q_DD = exp(log_q_DD);

  //--DECLARING DERIVED VALUES
  Type BMSY_DD = MSY_DD/UMSY_DD;
  Type SS_DD = So_DD * (1 - UMSY_DD);
  Type Spr_DD = (SS_DD * Alpha_DD/(1 - SS_DD) + wa_DD)/(1 - Rho_DD * SS_DD);
  Type DsprDu_DD = (Alpha_DD + Spr_DD * (1 + Rho_DD - 2 * Rho_DD * SS_DD))/((1 - Rho_DD * SS_DD) * (1 - SS_DD));
  DsprDu_DD += Alpha_DD * SS_DD/((1 - Rho_DD * SS_DD) * pow(1 - SS_DD, 2));
  DsprDu_DD -= Spr_DD/(1 - SS_DD);
  DsprDu_DD *= -So_DD;
  Type Arec_DD = 1/(pow(1 - UMSY_DD,2) * (Spr_DD + UMSY_DD * DsprDu_DD));
  Type Brec_DD = UMSY_DD * (Arec_DD * Spr_DD - 1/(1 - UMSY_DD))/MSY_DD;
  Type Spr0_DD = (So_DD * Alpha_DD/(1 - So_DD) + wa_DD)/(1 - Rho_DD * So_DD);
  Type Ro_DD = (Arec_DD * Spr0_DD - 1)/(Brec_DD * Spr0_DD);
  Type Bo_DD = Ro_DD * Spr0_DD;
  Type No_DD = Ro_DD/(1 - So_DD);

  //--DECLARING STORAGE VECTORS
  int ny_DDp = ny_DD + 1;
  int ny_DDk = ny_DD + k_DD;
  vector<Type> B_DD(ny_DDp);
  vector<Type> relB_DD(ny_DDp);
  vector<Type> N_DD(ny_DDp);
  vector<Type> R_DD(ny_DDk);

  vector<Type> Surv_DD(ny_DD);
  vector<Type> Cpred_DD(ny_DD);
  vector<Type> Sp_DD(ny_DD);
  vector<Type> U_DD(ny_DD);
  vector<Type> relU_DD(ny_DD);

  //--INITIALIZE
  B_DD(0) = Bo_DD;
  relB_DD(0) = B_DD(0)/BMSY_DD;
  N_DD(0) = No_DD;
  for(int tt=0;tt<k_DD;tt++) R_DD(tt) = Ro_DD;

  for(int tt=0; tt<ny_DD; tt++){
    U_DD(tt) = 1 - exp(-q_DD * E_hist(tt));
    relU_DD(tt) = U_DD(tt)/UMSY_DD;
    Surv_DD(tt) = So_DD * (1 - U_DD(tt));
    Cpred_DD(tt) = CppAD::CondExpGt(U_DD(tt) * B_DD(tt), Type(1e-15), U_DD(tt) * B_DD(tt), Type(1e-15));
    Sp_DD(tt) = B_DD(tt) - Cpred_DD(tt);

    R_DD(tt + k_DD) = Arec_DD * Sp_DD(tt)/(1 + Brec_DD * Sp_DD(tt));
    B_DD(tt + 1) = Surv_DD(tt) * (Alpha_DD * N_DD(tt) + Rho_DD * B_DD(tt)) + wa_DD * R_DD(tt + 1);
	relB_DD(tt + 1) = B_DD(tt + 1)/BMSY_DD;
    N_DD(tt + 1) = Surv_DD(tt) * N_DD(tt) + R_DD(tt + 1);
  }

  //--ARGUMENTS FOR NLL
  Type sigma_DD = calc_sigma(C_hist, Cpred_DD);

  // The following conditions must be met for positive values
  // of Arec_DD and Brec_DD, respectively:
  // umsy * DsprDu + Spr_DD > 0 and Arec_DD * Spr_DD * (1 - UMSY_DD) - 1 > 0
  // Thus, create a likelihood penalty of 1000 when either condition is not met
  Type penalty = CppAD::CondExpGt(Spr_DD + UMSY_DD * DsprDu_DD, Type(0), Type(0), Type(UMSY_DD * 1e3));
  penalty += CppAD::CondExpGt(Arec_DD * Spr_DD * (1 - UMSY_DD) - 1, Type(0), Type(0), Type(UMSY_DD * 1e3));

  // Objective function
  //creates storage for jnll and sets value to 0
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();

  for(int tt=0; tt<ny_DD; tt++){
    jnll_comp(0) -= dnorm(log(C_hist(tt)), log(Cpred_DD(tt)), sigma_DD, true);
  }
  //jnll_comp(1) -= dbeta(UMSY_DD, UMSYprior(0), UMSYprior(1), true);

  //Summing individual jnll and penalties
  Type jnll = jnll_comp.sum() + penalty;

  //-------REPORTING-------//
  Type TAC = UMSY_DD * B_DD(ny_DD);
  Type h = Arec_DD * Spr0_DD / (4 + Arec_DD * Spr0_DD);

  ADREPORT( UMSY_DD );
  ADREPORT( MSY_DD );
  ADREPORT( q_DD );
  ADREPORT( sigma_DD );
  REPORT( sigma_DD );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Arec_DD );
  REPORT( Brec_DD );
  REPORT( Spr_DD );
  REPORT( Spr0_DD );
  REPORT( DsprDu_DD );
  REPORT( Cpred_DD );
  REPORT( B_DD );
  REPORT( N_DD );
  REPORT( R_DD );
  REPORT( U_DD );
  REPORT( relU_DD );
  REPORT( relB_DD );
  REPORT( TAC );
  REPORT( h );
  REPORT( BMSY_DD );
  REPORT( Ro_DD );
  REPORT( No_DD );
  REPORT( Bo_DD );
  REPORT( penalty );

  return jnll;
//}
