#' Delay - Difference Stock Assessment in TMB with UMSY and MSY as leading parameters
#'
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded in TMB. The model
#' is conditioned on effort and estimates predicted catch.
#'
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param report Indicates whether report will be produced for Data object.
#' @return A numeric vector of TAC recommendations. If \code{report = TRUE}, a list of
#' model and TAC output is returned.
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on assumptions such as stationary productivity and
#' proportionality between the abundance index and real abundance.
#' Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers & Z. Siders. Zach Siders coded the TMB function
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses).
#' @export DD_TMB
#' @seealso \code{\link{DD_SS}}
#' @import TMB
#' @importFrom stats nlminb
#' @importFrom mvtnorm rmvnorm
#' @useDynLib MSE
DD_TMB <- function(x, Data, reps = 100, report = FALSE) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind"
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x,   ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYpriorpar <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3) # Prior for UMSY is that corresponding to F = 0.5 M with CV = 0.3
  UMSYprior <- c(alphaconv(UMSYpriorpar[1], prod(UMSYpriorpar)), betaconv(UMSYpriorpar[1], prod(UMSYpriorpar))) # Convert to beta parameters
  AvC <- mean(C_hist, na.rm = TRUE)
  data <- list(model = "DD", So_DD = So_DD, Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD,
               wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist, UMSYprior = UMSYprior)
  params <- list(logit_UMSY_DD = log(UMSYpriorpar[1]/(1 - UMSYpriorpar[1])),
                 log_MSY_DD = log(3 * AvC), log_q_DD = log(Data@Mort[x]))
  info <- list(data = data, params = params)

  obj <- MakeADFun(data = info$data, parameters = info$params, DLL = "MSE", silent = TRUE)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
  if(reps == 1) TAC <- obj$report()$TAC
  if(reps > 1) {
    SD <- sdreport(obj, getReportCovariance = FALSE)
    samps <- rmvnorm(reps, opt$par, round(SD$cov.fixed, 4))
    TAC <- rep(NA, reps)
    for (i in 1:reps) {
      params.new <- list(logit_UMSY_DD = samps[i, 1], log_MSY_DD = samps[i, 2],
                         log_q_DD = samps[i, 3])
      obj.samp <- MakeADFun(data = info$data, parameters = params.new, DLL = "MSE")
      TAC[i] <- obj.samp$report()$TAC
    }
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)

  if(report) {
    return(list(Rec = Rec, TAC = TAC, info = info, obj = obj, opt = opt,
                Data = Data, dependencies = dependencies))
  } else {
    return(Rec)
  }
}
class(DD_TMB) <- "MP"
