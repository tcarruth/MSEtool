#' State-Space Delay-Difference Stock Assessment in TMB with UMSY and MSY as leading parameters
#'
#' A delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded in TMB.
#' Recruitment deviations from the stock-recruit relationship are estimated as a
#' state space variable. The model is conditioned on effort and estimates
#' predicted catch.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data A data-limited methods data object.
#' @param ... Additional arguments (not currently used).
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
DD_SS <- function(x, Data, ...) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@CV_Cat, Data@Ind"
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
  ny <- length(C_hist)
  k <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)
  k[k > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho <- (wa[k + 2] - Winf)/(wa[k + 1] - Winf)
  Alpha <- Winf * (1 - Rho)
  S0 <- exp(-Data@Mort[x])  # get So survival rate
  wk <- wa[k]
  UMSYstart <- 1 - exp(-Data@Mort[x] * 0.5)
  AvC <- mean(C_hist, na.rm = TRUE)
  sigmaC <- max(0.05, sdconv(AvC, AvC * Data@CV_Cat[x]))
  data <- list(model = "DD_SS", S0 = S0, Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, E_hist = E_hist, C_hist = C_hist)
  params <- list(logit_UMSY = log(UMSYstart/(1 - UMSYstart)),
                 log_MSY = log(3 * AvC), log_q = log(Data@Mort[x]),
                 log_sigma = log(sigmaC), log_tau = log(0.3), log_rec_dev = rep(0, ny - k))
  info <- list(Year = Year, data = data, params = params, sigma = sigmaC, I_hist = I_hist)

  obj <- MakeADFun(data = info$data, parameters = info$params, random = "log_rec_dev",
                   map = list(log_sigma = factor(NA)), DLL = "MSEtool", silent = TRUE)
  opt <- optimize_TMB_model(obj)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || opt$convergence != 0 || is.character(SD)) {
    Assessment <- new("Assessment", Model = "DD_SS",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
    warning("Model did not properly converge. Check TMB objects (slots names: opt and SD).")
  } else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearplusk <- c(Year, max(Year) + 1:k)
    Yearrandom <- seq(Year[1] + k, max(Year))

    Assessment <- new("Assessment", Model = "DD_SS",
                      UMSY = report$UMSY, MSY = report$MSY, BMSY = report$BMSY,
                      B0 = report$B0, R0 = report$R0, N0 = report$N0,
                      SSB0 = report$B0, h = report$h,
                      U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$B0, names = Yearplusone),
                      SSB = structure(report$B, names = Yearplusone),
                      SSB_SSBMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      SSB_SSB0 = structure(report$B/report$B0, names = Yearplusone),
                      R = structure(report$R, names = Yearplusk),
                      N = structure(report$N, names = Yearplusone),
                      Obs_Catch = structure(C_hist, names = Year),
                      Catch = structure(report$Cpred, names = Year),
                      Random = structure(SD$par.random, names = Yearrandom),
                      Random_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$jnll_comp), names = c("Total", "Catch", "Random")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[6],
                      SE_B_BMSY_final = SD$sd[7], SE_B_B0_final = SD$sd[8],
                      SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(DD_SS) <- "Assess"

