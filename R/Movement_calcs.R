
# Movement matrix calculation functions
#' Calculates movement matrices from user inputs
#'
#' @description A function for calculating a movement matrix from user specified unfished stock biomass fraction in each area
#' @param OM Operating model, an object of class \linkS4class{OM}.
#' @param dist A vector of length \code{nareas} of fractions of unfished stock in each area
#' @param prob Mean probability of staying across all areas (single value) or a vector of the probability of individuals staying in each area (same length as dist)
#' @param distE Logit (normal) St.Dev error for sampling stock fractions from the fracs vector
#' @param probE Logit (normal) St.Dev error for sampling desired probability of staying either by area (prob is same length as dist) or the mean probability of staying (prob is a single number)
#' @param prob2 Optional vector as long as prob and dist. Upper bounds on uniform sampling of probability of staying, lower bound is prob.
#' @author T. Carruthers
#' @export simmov
#' @import TMB
#' @useDynLib MSEtool
#' @examples
#' movOM_5areas<-simmov(testOM,dist=c(0.01,0.1,0.2,0.3,0.39),prob=c(0.1,0.6,0.6,0.7,0.9))
#' movOM_5areas@cpars$mov[1,1,,] # sim 1, age 1, movement from areas to areas
#'
simmov<-function(OM,dist=c(0.1,0.2,0.3,0.4),prob=0.5,distE=0.1,probE=0.1,prob2=NA){

  logit<-function(p)log(p/(1-p))
  ilogit<-function(x)exp(x)/(1+exp(x))
  ilogitm<-function(x)exp(x)/apply(exp(x),1,sum)

  nareas<-length(dist)
  nsim<-OM@nsim
  maxage<-OM@maxage

  mov<-array(NA,c(nsim,maxage,nareas,nareas))

  dist_s<-ilogitm(t(matrix(rnorm(nareas*nsim,log(dist),distE),nrow=nareas)))

  if(length(prob)==1){

    prob_s<-ilogit(matrix(rnorm(nsim,logit(prob),probE),nrow=nsim))

  }else if(length(prob)==length(dist)&is.na(prob2)){

    prob_s<-ilogit(t(matrix(rnorm(nareas*nsim,logit(prob),probE),nrow=nareas)))

  }else if(length(prob)==length(dist)&length(prob)==length(prob2)){

    prob_s<-t(matrix(runif(nareas*nsim,prob,prob2),nrow=nareas))

  }else{

    stop("Error: either prob wasn't of length 1, or prob wasn't of length dist or prob 2 wasn't the same length as prob and dist.
            You have three options:
            (1) provide one value for prob which represents mean probability of staying across all areas sampled for each simulation with probE logit error
            (2) provide nareas values of prob which represent probability of staying across all areas sampled for each simulation with probE logit error
            (3) provide nareas values of prob and prob2 which are the upper and lower bounds for sampling uniform probability of staying for each area")
  }

  for(i in 1:nsim){

    movt<-makemov(fracs=dist_s[i,],prob=prob_s[i,])
    mov[i,,,]<-rep(movt,each=maxage)
  } # nsim

  OM@cpars$mov<-mov
  OM

}
# myOM<-simmov(testOM, dist=c(0.1,0.8,0.05,0.05),prob=c(0.6,0.5,0.7,0.85))


#' Calculates movement matrices from user inputs for fraction in each area (fracs) and probability of staying in areas (prob)
#'
#' @description A function for calculating a movement matrix from user specified unfished stock biomass fraction in each area
#' @param fracs A vector nareas long of fractions of unfished stock biomass in each area
#' @param prob A vector of the probability of individuals staying in each area or a single value for the mean probability of staying among all areas
#' @author T. Carruthers
#' @export makemov
#' @import TMB
#' @importFrom stats nlminb
#' @useDynLib MSEtool
#' @examples
#' makemov(fracs=c(0.1,0.5,0.2,0.2),prob=c(0.9,0.5,0.3,0.8))
makemov<-function(fracs=c(0.1,0.2,0.3,0.4),prob=c(0.5,0.8,0.9,0.95)){

  nareas<-length(fracs)

  if(length(prob)==1) data <- list(model = "grav",fracs = fracs, prob = prob, nareas = nareas)
  if(length(prob)==nareas) data <- list(model = "grav_Pbyarea",fracs = fracs, prob = prob, nareas = nareas)

  if(length(prob)==1)params <- list(log_visc = 0,log_grav = rep(0,nareas-1))
  if(length(prob)==nareas)params <- list(log_visc = rep(0,nareas),log_grav = rep(0,nareas-1))

  info <- list(data = data, params = params)

  obj <- MakeADFun(data = info$data, parameters = info$params, DLL = "MSEtool", silent = TRUE)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)

  if(length(prob)==1)params.new<-list(log_visc=opt$par[1],log_grav=opt$par[2:nareas])
  if(length(prob)==nareas)params.new<-list(log_visc=opt$par[1:nareas],log_grav=opt$par[nareas+(1:(nareas-1))])

  info <- list(data = data, params = params.new)
  obj.new<-MakeADFun(data = info$data, parameters = info$params, DLL = "MSEtool", silent = TRUE)

  obj.new$report()$mov
  #validateTMB(obj.new)

}

#' Checks the TMB equations again estimated movement matrices and also checks fit
#'
#' @description A function for calculating a movement matrix from user specified unfished stock biomass fraction in each area
#' @param obj A list object arising from MakeADFun from either grav.h or grav_Pbyarea.h
#' @author T. Carruthers
#' @export validateTMB
#' @import TMB
validateTMB<-function(obj){

  log_grav<-obj$report()$log_grav
  log_visc<-obj$report()$log_visc
  nareas<-length(log_grav)+1

  grav<-array(0,c(nareas,nareas))
  grav[,2:nareas]<-rep(log_grav,each=nareas)
  grav[cbind(1:nareas,1:nareas)]<-grav[cbind(1:nareas,1:nareas)]+log_visc

  mov<-exp(grav)/apply(exp(grav),1,sum)

  idist<-rep(1/nareas,nareas)
  for(i in 1:50)idist<-apply(idist*mov,2,sum)

  print(obj$report()$idist)
  print(idist)
  print(obj$report()$fracs)

  print(obj$report()$mov)
  print(mov)

}



