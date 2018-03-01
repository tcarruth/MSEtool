
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace density;

  DATA_VECTOR( fracs );
  DATA_SCALAR( muprob );
  DATA_INTEGER( nareas );

  PARAMETER_VECTOR( log_grav );
  PARAMETER( log_visc );

  vector<Type> log_grav_p = log_grav; // this is a vector, not sure what will happen here
  Type log_visc_p = log_visc;
  Type gsum;
  Type psum;
  Type sigma;
  Type obj;

  // -- Declarations
  matrix<Type> grav(nareas,nareas);
  matrix<Type> mov(nareas,nareas);
  matrix<Type> transN(nareas,nareas);
  vector<Type> idist(nareas);

  // Map out gravity terms accounting for viscosity, gravity of area 1 is fixed to zero
  for(int af=1; af<nareas; af++){ // area from
    grav(af,1) = 0;
    for(int at=2; af<nareas; af++){ // area to
      grav(af,at) = log_grav_p(at-1);
    }
    grav(af,af)+=grav(af,af)+log_visc_p; // add viscosity
  }

  // Calculate logit fractions (movement to area from area)
  for(int af=1; af<nareas; af++){
    gsum=0.0;
    for(int at=1; af<nareas; af++){
      gsum+=exp(grav(af,at)); // sum up gravity terms
    }
    for(int at=1; af<nareas; af++){
      mov(af,at)+=exp(grav(af,at))/gsum; // calculate logit mov probs by row (area from)
    }
  }

  // Run a convergence to a stable distribution
  for(int af=1; af<nareas; af++){
    idist(af)=1.0/nareas;
  }

  for(int tt =1; tt<20; tt++){
    for(int af=1; af<nareas; af++){
      for(int at=1; at<nareas; at++){
        transN(af,at)=idist(af)*mov(af,at);
      }
    }
    for(int at=1; at<nareas; at++){
      gsum=0.0;
      for(int af=1; af<nareas; af++){
        gsum+=transN(af,at);
      }
      idist(at)=gsum;
    }
  }

  psum=0.0;
  for(int aa=1; aa<nareas; aa++){
    psum+=mov(aa,aa)/nareas;
  }

  sigma = 0.1000;

  gsum=0.0;
  for(int aa=1; aa<nareas; aa++){
    gsum-=dnorm(log(idist(aa)), log(fracs(aa)), sigma, true);
  }

  gsum-=dnorm(log(psum), log(muprob), sigma, true);

  obj = gsum;

  //-------REPORTING-------//
  ADREPORT( log_grav_p );
  ADREPORT( log_visc_p );
  REPORT( idist );
  REPORT( transN );

  return obj;
//}
