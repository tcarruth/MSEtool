template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
    //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
    Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
    if(give_log) return logres; else return exp(logres);
}

