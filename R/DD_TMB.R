#' Delay - Difference Stock Assessment in TMB with UMSY and MSY as leading parameters
#'
#' A simple delay-difference assessment that estimates the TAC using a
#' time-series of catches and a relative abundance index and coded in TMB. The model
#' is conditioned on effort and estimates predicted catch.
#'
#' @param x An index for the objects in \code{Data} when running in closed loop simulation.
#' Otherwise, equals to 1 When running an assessment interactively.
#' @param Data An object of class \linkS4class{Data}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing objects and output
#' from TMB.
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on assumptions such as stationary productivity and
#' proportionality between the abundance index and real abundance.
#' Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers & Z. Siders. Zach Siders coded the TMB function
#' @references
#' Carruthers, T, Walters, C.J,, and McAllister, M.K. 2012. Evaluating methods that classify
#' fisheries stock status using only fisheries catch data. Fisheries Research 119-120:66-79.
#'
#' Hilborn, R., and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York.
#' @export DD_TMB
#' @seealso \code{\link{DD_SS}}
#' @import TMB
#' @importFrom stats nlminb
#' @useDynLib MSEtool
DD_TMB <- function(x, Data, ...) {
  dependencies = "Data@vbLinf, Data@vbK, Data@vbt0, Data@Mort, Data@wla, Data@wlb, Data@Cat, Data@Ind, Data@L50"
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
  UMSY_start <- 1 - exp(-Data@Mort[x] * 0.5)
  AvC <- mean(C_hist, na.rm = TRUE)
  data <- list(model = "DD", S0 = S0, Alpha = Alpha, Rho = Rho, ny = ny, k = k,
               wk = wk, E_hist = E_hist, C_hist = C_hist)
  params <- list(logit_UMSY = log(UMSY_start/(1 - UMSY_start)),
                 log_MSY = log(3 * AvC), log_q = log(Data@Mort[x]))
  info <- list(Year = Year, data = data, params = params, I_hist = I_hist)

  obj <- MakeADFun(data = info$data, parameters = info$params, DLL = "MSEtool", silent = TRUE)
  opt <- optimize_TMB_model(obj)
  SD <- get_sdreport(obj, opt)
  report <- obj$report(obj$env$last.par.best)

  if(is.character(opt) || opt$convergence != 0 || is.character(SD)) {
    Assessment <- new("Assessment", Model = "DD_TMB",
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
    warning("Model did not converge. Check TMB objects (slots names: obj, opt, and SD).")
  } else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearplusk <- c(Year, max(Year) + 1:k)

    Assessment <- new("Assessment", Model = "DD_TMB",
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
                      NLL = structure(c(opt$objective, report$nll), names = c("Total", "Catch")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2],
                      SE_U_UMSY_final = SD$sd[5], SE_B_BMSY_final = SD$sd[6],
                      SE_B_B0_final = SD$sd[7], info = info, obj = obj, opt = opt,
                      SD = SD, TMB_report = report, dependencies = dependencies, Data = Data)
  }

  return(Assessment)
}
class(DD_TMB) <- "Assess"
