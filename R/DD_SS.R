#' State-Space Delay-Difference Stock Assessment in TMB with UMSY and MSY as leading parameters
#'
#' A delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded in TMB.
#' Recruitment deviations from the stock-recruit relationship are estimated as a
#' state space variable. The model is conditioned on effort and estimates
#' predicted catch.
#'
#' @param Data A data-limited methods data object
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @note Similar to many other assessment
#' models it depends on assumptions such as stationary productivity and
#' proportionality between the abundance index and real abundance.
#' Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' The observation error standard deviation is fixed according to the
#' presumed CV of the catch and the process error (recruitment) standard deviation
#' is estimated.
#' @author Q. Huynh
#' @references
#' Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating methods that classify
#' fisheries stock status using only fisheries catch data. Fisheries Research 119-120:66-79.
#'
#' Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York.
#' @export DD_SS
#' @seealso \code{\link{DD_TMB}}
#' @import TMB
#' @importFrom stats nlminb
#' @useDynLib MSEtool
DD_SS <- function(Data) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@CV_Cat, Data@Ind"
  x <- 1 # Legacy of Data structure
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x],  Data@L50[x])
  a50V <- max(a50V, 1)
  ystart <- which(!is.na(Data@Cat[x, ] + Data@Ind[x,   ]))[1]
  yind <- ystart:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  I_hist <- Data@Ind[x, yind]
  E_hist <- C_hist/I_hist
  if(any(is.na(E_hist))) stop("Missing values in catch and index in Data object.")
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYpriorpar <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3) # Prior for UMSY is that corresponding to F = 0.5 M with CV = 0.3
  #UMSYprior <- c(alphaconv(UMSYpriorpar[1], prod(UMSYpriorpar)), betaconv(UMSYpriorpar[1], prod(UMSYpriorpar))) # Convert to beta parameters
  AvC <- mean(C_hist, na.rm = TRUE)
  sigmaC <- max(0.05, sdconv(AvC, AvC * Data@CV_Cat[x]))
  data <- list(model = "DD_SS", So_DD = So_DD, Alpha_DD = Alpha_DD, Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD,
               wa_DD = wa_DD, E_hist = E_hist, C_hist = C_hist)
  params <- list(logit_UMSY_DD = log(UMSYpriorpar[1]/(1 - UMSYpriorpar[1])),
                 log_MSY_DD = log(3 * AvC), log_q_DD = log(Data@Mort[x]),
                 log_sigma_DD = log(sigmaC),
                 log_tau_DD = log(0.3), log_rec_dev = rep(0, ny_DD - k_DD))
  info <- list(Year = Year, data = data, params = params, sigma = sigmaC, I_hist = I_hist)

  obj <- MakeADFun(data = info$data, parameters = info$params, random = "log_rec_dev",
                   map = list(log_sigma_DD = factor(NA)), DLL = "MSEtool", silent = TRUE)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
  SD <- sdreport(obj)

  Assessment <- return_Assessment()
  return(Assessment)
}
class(DD_SS) <- "Assess"

