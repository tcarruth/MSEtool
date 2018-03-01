
# Movement matrix calculation functions


#' Calculates movement matrices from user inputs
#'
#' @description A function for calculating a movement matrix from user specified unfished stock biomass fraction in each area
#' @param fracs A vector nareas long of fractions of unfished stock biomass in each area
#' @param visc A vector of the probability of individuals staying in each area
#' @param efracs A log error term for sampling stock fractions from the fracs vector
#' @param evisc A vector of log error terms for sampling visc by area
#' @author T. Carruthers
#' @export makemov
#' @import TMB
#' @importFrom stats nlminb
#' @importFrom mvtnorm rmvnorm
#' @useDynLib MSEtool
makemov<-function(fracs=rep(0.2,5),muprob=0.5,efracs=0.1,evisc=rep(0.1,5),OM=NULL){

  nareas<-length(fracs)
  data <- list(model = "grav",fracs = fracs, muprob = muprob, nareas = nareas)

  params <- list(log_visc = 0,log_grav = rep(0,nareas-1))
  info <- list(data = data, params = params)

  obj <- MakeADFun(data = info$data, parameters = info$params, DLL = "MSEtool", silent = TRUE)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)

  if(reps == 1) TAC <- obj$report()$TAC
  if(reps > 1) {
    SD <- sdreport(obj, getReportCovariance = FALSE)
    samps <- rmvnorm(reps, opt$par, round(SD$cov.fixed, 4))
    TAC <- rep(NA, reps)
    for (i in 1:reps) {
      params.new <- list(logit_UMSY_DD = samps[i, 1], log_MSY_DD = samps[i, 2],
                         log_q_DD = samps[i, 3])
      obj.samp <- MakeADFun(data = info$data, parameters = params.new, DLL = "MSEtool")
      TAC[i] <- obj.samp$report()$TAC
    }
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)

}


