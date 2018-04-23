#' Statistical catch-at-age (SCA) model with UMSY and MSY as leading parameters
#'
#' A generic statistical catch-at-age model. The parameterization by Forrest et al. (2008)
#' is used which generally results in faster convergence.
#'
#' @param x A position in the Data object (by default, equal to one for
#' assessments).
#' @param Data An object of class Data.
#' @param U_begin Indicates how to handle the equilibrium harvest rate prior to the first
#' year of model. \code{"virgin"} fixes \code{U = 0} prior to year 1 while \code{"est"} will
#' estimate it.
#' @param vulnerability Whether estimated vulnerability is logistic or dome (double-normal).
#' See details for parameterization.
#' @param SR Stock-recruit function (currently only supports Beverton-Holt).
#' @param CAA_multiplier Numeric for data weighting of catch-at-age matrix. See details.
#' @param start Optional list of starting values. See details.
#' @param silent (TRUE/FALSE) Passed to \code{\link[TMB]{MakeADFun}}, whether TMB
#' will print trace information during optimization. Used for dignostics for model convergence.
#' @param control A named list of parameters regarding optimization to be passed to
#' \code{\link[stats]{nlminb}}.
#' @param ... Other arguments to be passed.
#' @return An object of class \linkS4class{Assessment}.
#' @note Convergence problems may occur if perfect data are used (i.e., \code{Data@@CV_Ind}
#' is very low).
#' @details
#' For the SCA model, the basic data inputs are catch (by weight), index
#' (by weight/biomass), and catch-at-age matrix (by numbers). Catches are
#' assumed to be known perfectly (the harvest rate in a given year is the ratio of the observed
#' catch to the vulnerable biomass at the beginning of the year), while the observation error
#' of the index is fixed based on the value of \code{CV_Ind} in the Data object.
#'
#' The annual sample sizes of the catch-at-age matrix is provided to the model (used in the
#' likelihood for catch-at-age, assuming a multinomial distribution),
#' and is manipulated via argument \code{CAA_multiplier}. This argument is
#' interpreted in two different ways depending on the value provided.
#' If \code{CAA_multiplier > 1}, then this value will cap the annual sample sizes
#' to that number. If \code{CAA_multiplier <= 1}, then all the annual samples sizes
#' will be re-scaled by that number. By default, sample sizes are capped at 50.
#'
#' Recruitment deviations are estimated beginning in the year when age composition are available.
#'
#' Vulnerability can be specified to be either logistic or dome. If logistic, then the model
#' vector \code{vul_par} is of length 2, containing the ages of 50% and 95% vulnerability,
#' respectively. The age of 95% vulnerability is an offset, i.e.,
#' \code{vul_95 = vul_par[1] + exp(vul_par[2])}.
#'
#' With dome vulnerability, a double normal parameterization is used, where \code{vul_par}
#' is an estimated vector of length 4:
#' \itemize{
#' \item \code{vulpar[1]} \code{log(sd_asc)}, where sd_asc of the normal distribution function for the ascending limb
#' \item \code{vulpar[2]} \code{mu_asc}, mean of the normal distribution function for the ascending limb
#' \item \code{vulpar[3]} \code{mu_des}, mean of the normal distribution function for the descending limb. To ensure
#' \code{mu_desc > mu_asc}, this is parameterized as an offset, i.e., \code{mu_desc = mu_asc + exp(vulpar[3])}.
#' \item \code{vulpar[4]} \code{log(sd_des)}, where sd_des of the normal distribution function for the descending limb.
#' }
#'
#' For \code{start}, a named list of starting values of estimates for:
#' \itemize{
#' \item \code{UMSY}
#' \item \code{MSY}
#' \item \code{U_equilibrium} (typically 0). See \code{U_begin} argument.
#' \item \code{vul_par} (length 2 vector for logistic or length 4 for dome).
#' \item \code{tau} the standard deviation of the log-recruitment deviations. The observation
#' error of the index is fixed based on the value of \code{CV_Ind} in the Data object, and
#' weighting of the catch-at-age is handled by argument \code{CAA_multiplier}.
#' }
#' @author Q. Huynh
#' @references
#' Forrest, R.E., Martell, S.J.D., Melnychuk, M.C., and Walters, C.J. 2008.
#' An age-structured model with leading management parameters, incorporating
#' age-specific selectivity and maturity. Canadian Journal of Fisheries and Aquatic
#' Science 65:286-296.
#' @export
SCA <- function(x = 1, Data, U_begin = c("virgin", "est"), vulnerability = c("logistic", "dome"),
                SR = c("BH", "Ricker"), CAA_multiplier = 50, start = NULL,
                silent = TRUE, control = list(), ...) {
  dependencies = ""
  U_begin <- match.arg(U_begin)
  vulnerability <- match.arg(vulnerability)
  SR <- match.arg(SR)
  yind <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- yind:length(Data@Cat[x, ])
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if(any(is.na(C_hist) | C_hist < 0)) warning("Error. Catch time series is not complete.")
  I_hist <- Data@Ind[x, yind]
  Data <- expand_comp_matrix(Data, "CAA") # Make sure dimensions of CAA match that in catch (nyears).
  CAA_hist <- Data@CAA[x, yind, ]

  CAA_n_nominal <- rowSums(CAA_hist)
  if(CAA_multiplier <= 1) {
    CAA_n_rescale <- CAA_multiplier * CAA_n_nominal
  } else CAA_n_rescale <- pmin(CAA_multiplier, CAA_n_nominal)

  n_y <- length(C_hist)
  max_age <- Data@MaxAge
  M <- rep(Data@Mort[x], max_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(1:max_age) - t0)))
  Wa <- a * La ^ b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50+0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - A50)))

  CAA_mode <- which.max(colSums(CAA_hist, na.rm = TRUE))
  if(vulnerability == "logistic") {
    vul_par <- c(CAA_mode-1, log(1)) # 50 and log(95%-offset) vulnerability respectively
  }
  if(vulnerability == "dome") {
    # double normal: logsd(ascending), mean(asc), mean(desc)-logoffset, sd(desc)
    vul_par <- c(log(1), CAA_mode, log(0.5), log(5))
  }

  UMSYstart <- 1 - exp(-Data@Mort[x] * 0.5)
  AvC <- mean(C_hist, na.rm = TRUE)
  sigma_I <- max(0.05, sdconv(1, Data@CV_Ind[x])) # can't fit to perfect data

  data <- list(model = "SCA", C_hist = C_hist, I_hist = I_hist,
               CAA_hist = t(apply(CAA_hist, 1, function(x) x/sum(x))),
               CAA_n = CAA_n_rescale, n_y = n_y, max_age = max_age, M = M,
               weight = Wa, mat = mat_age, vul_type = vulnerability,
               est_rec_dev = as.integer(random_map(CAA_n_nominal)))

  params <- list(logit_UMSY = log(UMSYstart[1]/(1 - UMSYstart[1])),
                 log_MSY = log(3 * AvC), U_equilibrium = 0, vul_par = vul_par,
                 log_sigma = log(sigma_I), log_tau = log(0.5),
                 log_rec_dev = rep(0, n_y - 1))
  info <- list(Year = Data@Year, data = data, params = params, control = control)

  map <- list(log_sigma = factor(NA))
  if(U_begin == "virgin") map$U_equilibrium <- factor(NA)
  if(any(is.na(CAA_n_nominal))) {
    map$log_rec_dev <- random_map(CAA_n_nominal)
  }
  obj <- MakeADFun(data = info$data, parameters = info$params,
                   map = map, random = "log_rec_dev",
                   DLL = "MSEtool", silent = silent)
  opt.time <- system.time(opt <- optimize_TMB_model(obj, control))[3]
  SD.time <- system.time(SD <- get_sdreport(obj, opt))[3]
  report <- obj$report(obj$env$last.par.best)
  if(is.character(opt) || is.character(SD)) {
    Assessment <- new("Assessment", Model = "SCA", info = info,
                      obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  else {
    Yearplusone <- c(Year, max(Year) + 1)
    Yearrandom <- seq(Year[2], max(Year))
    Assessment <- new("Assessment", Model = "SCA", UMSY = report$UMSY,
                      MSY = report$MSY, BMSY = report$BMSY, SSBMSY = report$EMSY,
                      VBMSY = report$VBMSY, B0 = report$B0, R0 = report$R0, N0 = report$N0,
                      SSB0 = report$E0, VB0 = report$VB0,
                      h = report$h, U = structure(report$U, names = Year),
                      U_UMSY = structure(report$U/report$UMSY, names = Year),
                      B = structure(report$B, names = Yearplusone),
                      B_BMSY = structure(report$B/report$BMSY, names = Yearplusone),
                      B_B0 = structure(report$B/report$B0, names = Yearplusone),
                      SSB = structure(report$E, names = Yearplusone),
                      SSB_SSBMSY = structure(report$E/report$EMSY, names = Yearplusone),
                      SSB_SSB0 = structure(report$E/report$E0, names = Yearplusone),
                      VB = structure(report$VB, names = Yearplusone),
                      VB_VBMSY = structure(report$VB/report$VBMSY, names = Yearplusone),
                      VB_VB0 = structure(report$VB/report$VB0, names = Yearplusone),
                      R = structure(report$R, names = Yearplusone),
                      N = structure(rowSums(report$N), names = Yearplusone),
                      N_at_age = report$N,
                      Selectivity = matrix(report$vul, nrow = length(Year),
                                           ncol = max_age, byrow = TRUE),
                      Obs_Catch = structure(C_hist, names = Year),
                      Obs_Index = structure(I_hist, names = Year),
                      Obs_C_at_age = CAA_hist,
                      Index = structure(report$Ipred, names = Year),
                      C_at_age = report$CAApred,
                      Random = structure(SD$par.random, names = Yearrandom),
                      Random_type = "log-Recruitment deviations",
                      NLL = structure(c(opt$objective, report$nll_comp),
                                      names = c("Total", "Index", "CAA", "Random")),
                      SE_UMSY = SD$sd[1], SE_MSY = SD$sd[2], SE_U_UMSY_final = SD$sd[6],
                      SE_B_BMSY_final = SD$sd[7], SE_B_B0_final = SD$sd[8],
                      SE_SSB_SSBMSY_final = SD$sd[9], SE_SSB_SSB0_final = SD$sd[10],
                      SE_VB_VBMSY_final = SD$sd[11], SE_VB_VB0_final = SD$sd[12],
                      SE_Random = structure(sqrt(SD$diag.cov.random), names = Yearrandom),
                      info = info, obj = obj, opt = opt, SD = SD, TMB_report = report,
                      dependencies = dependencies, Data = Data)
  }
  return(Assessment)
}
class(SCA) <- "Assess"

