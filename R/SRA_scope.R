


#' Stock-reduction analysis (SRA) for conditioning operating models
#'
#' Intended for conditioning operating models for data-limited stocks. From a historical time series of total catch or effort, and potentially
#' age/length compositions and multiple indices of abundance, the SRA returns a range of values for depletion, selectivity,
#' unfished recruitment (R0), historical fishing effort, and recruitment deviations for the operating model. This is done by sampling life history parameters
#' provided by the user and fitting to the data in a statistical catch-at-age model (with the predicted catch equal to the observed catch).
#' This function is intended to generate a range of potential depletion scenarios that could be supported from sparse data.
#' Either a full catch (conditioned on catch) or effort (conditioned on effort) time series is needed but missing data (as NAs) are allowed for all other data types.
#'
#' @param OM An object of class \linkS4class{OM} that specifies natural mortality (M), growth (Linf, K, t0, a, b), stock-recruitment relationship,
#' steepness, maturity parameters (L50 and L50_95), standard deviation of recruitment variability (Perr), as well as index uncertainty (Iobs).
#' @param data A list of data inputs. See Data section below.
#' @param condition String to indicate whether the SRA model is conditioned on catch or effort.
#' @param selectivity A character vector of length nfleet to indicate \code{"logistic"} or \code{"dome"} selectivity for each fleet in \code{Chist}.
#' @param s_selectivity A vector of length nsurvey to indicate \code{"logistic"} or \code{"dome"} selectivity for each survey in \code{Index}. Use a number
#' for an age-specific index.
#' Only used if any of the corresponding entries of \code{data$I_type = "est"} or if a number is specified here.
#' @param LWT A named list of likelihood weights for the SRA model. See details.
#' @param comp_like A string indicating either "multinomial" (default) or "lognormal" distributions for the composition data.
#' @param ESS If \code{comp_like = "multinomial"}, a numeric vector of length two to cap the maximum effective samples size of the age and length compositions,
#' respectively, for the multinomial likelihood function. The effective sample size of an age or length composition sample is the minimum of ESS or the number of observations
#' (sum across columns). For more flexibility, set ESS to be very large and alter the arrays as needed.
#' @param max_F The maximum F for any fleet in the scoping model (higher F's in the model are penalized in the objective function).
#' @param cores Integer for the number of CPU cores for the stock reduction analysis.
#' @param integrate Logical, whether to treat recruitment deviations as penalized parameters (FALSE) or random effects (TRUE).
#' @param mean_fit Logical, whether to run an additional with mean values of life history parameters from the OM.
#' @param sims A logical vector of length \code{OM@@nsim} or a numberic vector indicating which simulations to keep.
#' @param drop_nonconv Logical, whether to drop non-converged fits of the SRA model.
#' @param ... Other arguments to pass in for starting values of parameters and fixing parameters. See details.
#' @return An object of class \linkS4class{SRA}, including the updated operating model object.
#'
#' @section Data:
#' One of indices, age compositions, or length compositions should be provided in addition to the historical catch or effort. Not all arguments
#' are needed to run the model (some have defaults, while others are ignored if not applicable depending on the data provided).
#'
#' The \code{data} list can include:
#'
#' \itemize{
#' \item Chist - A vector of historical catch, should be of length OM@@nyears. If there are multiple fleets: a matrix of OM@@nyears rows and nfleet columns.
#' Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).
#' \item Ehist - A vector of historical effort, should be of length OM@@nyears (see also \code{E_eq}).
#' \item Index - A vector of values of an index (of length OM@@nyears). If there are multiple surveys: a matrix of historical indices of abundances, with rows
#' indexing years and columns indexing surveys. Age-specific indices should be numbers-specific while all others are weight-based.
#' \item I_sd - A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
#' If not provided, this function will use values from \code{OM@@Iobs}.
#' \item I_type - A character vector of length nsurvey to indicate the type of biomass for which each index follows. Either \code{"B"} for
#' total biomass, or \code{"SSB"} for spawning biomass. If not provided, "B" is used. Use numbers if the index corresponds to a fleet in \code{Chist}.
#' Use \code{"est"} to set survey selectivity to be an independent component of the model, i.e., as an age-specific index or estimated separatel.
#' Note, this generally requires age \code{s_CAA} or length \code{s_CAL}compositions.
#' \item CAA - Fishery age composition matrix with nyears rows and OM@@maxage columns. If multiple fleets: an array with dimension: nyears, OM@@maxage, and nfleets.
#' \item CAL - Fishery Length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets: an array with dimension: nyears,
#' length bins, and nfleets.
#' \item ML - A vector of fishery mean length observations (length OM@@nyears), or if multiple fleets: matrix of dimension: nyears and nfleets. Generally, should not
#' be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.
#' \item ML_sd - The standard deviation (normal distribution) of the observed mean lengths. If there are multiple fleets, a vector of length nfleet.
#' If not provided, default value is \code{0.1 * mean(ML)}.
#' \item s_CAA - Survey age composition data, an array of dimension nyears, maxage, nsurvey.
#' \item s_CAL - Survey length composition data, an array of dimension nyears, length(length_bin), nsurvey.
#' \item length_bin - A vector for the midpoints of the length bins for \code{CAL} and \code{s_CAL}. All bin widths should be equal in size.
#' \item C_eq - A numeric vector of length nfleet for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
#' Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.
#' \item E_eq - The equilibrium effort for each fleet in \code{Ehist} prior to the first year of the operating model.
#' Zero (default) implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.
#' \item abs_I - Optional, an integer vector to indicate which indices are in absolute magnitude. Use 1 to set q = 1, otherwise use 0 to estimate q.
#' \item I_basis - Optional, an integer vector to indicate whether indices are biomass based (1) or abundance-based (0). By default, all are biomass-based.
#' }
#'
#' Selectivity is fixed to values sampled from \code{OM} if no age or length compositions are provided.
#'
#' @details
#' For \code{SRA_scope}, additional arguments can be passed to the model via \code{...}:
#'
#' \itemize{
#' \item vul_par: A matrix of 3 rows and nfleet columns for starting values for fleet selectivity. The three rows correspond
#' to L5 (length of 5 percent selectivity), LFS (length of full selectivity) and Vmaxlen (selectivity at length Linf). By default,
#' the starting values are values from the OM object.
#' \item s_vul_par: A matrix of 3 rows and nsurvey columns for starting values for fleet selectivity. Same setup as vul_par.
#' \item map_vul_par: The map argument for vul_par in TMB, see \link[TMB]{MakeADFun}, which indicates whether selectivity parameters are fixed
#' or estimated. A matrix of the same dimension as vul_par. If an entry is \code{NA}, the corresponding parameter is fixed in the model to the starting
#' value. Otherwise, an integer for each independent parameter. By default, selectivity is fixed if there are no age or length composition for that fleet
#' or survey, otherwise estimated.
#' \item map_s_vul_par: The map argument for the survey selectivity parameters (same dimension as s_vul_par).
#' \item map_log_rec_dev: A vector of length OM@@nyears that indexes which recruitment deviates are fixed (using NA) or estimated (a separate integer).
#' }
#'
#' Survey selectivity is estimated only if \code{s_CAA} or \code{s_CAL} is provided. Otherwise, the selectivity should
#' be mirrored to a fleet (vulnerable biomass selectivity) or indexed to total or spawning biomass (see \code{I_type}).
#'
#' \code{LWT} is an optional named list containing the likelihood weights (values > 0) with the possible options:
#' \itemize{
#' \item Chist: A vector of length nfleet.
#' \item Index: A vector of length nsurvey.
#' \item CAA, CAL, ML, C_eq: A vector of length nfleet for each.
#' \item s_CAA, s_CAL: A vector of length nsurvey for each.
#' }
#'
#' By default, all likelihood weights are equal to one if not specified by the user. Weighting for CAA and CAL can also be adjusted by changing the
#' multinomial sample size. For \code{CAA}, \code{CAL}, \code{s_CAA}, and \code{s_CAL}, the arrays should be set up so that
#' the annual number of observations (summed over columns) should be equal to the presumed effective sample size. Argument \code{ESS} provides a shortcut
#' to cap the the effective sample size.
#'
#' \code{plot_SRA_scope} is now deprecated in favor of \link{plot.SRA}.
#'
#' Parameters that were used in the fitting model are placed in objects in \code{OM@@cpars}.
#'
#' \code{Sub_cpars} is a convenient function to subset simulations
#' for the operating model, for example, to remove simulations from unconverged model fits or outlier simulations.
#'
#' @note If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the SRA as well.
#' Time-varying life history parameters create ambiguity in the calculation and interpretation of depletion and reference points in \link[DLMtool]{runMSE}.
#' See section D.5 of \code{DLMtool::userguide()}.
#'
#' The easiest way to turn off time-varying growth/M is by setting: \code{OM@@Msd <- OM@@Linfsd <- OM@@Ksd <- c(0, 0)}.
#'
#' @author Q. Huynh
#' @seealso \link{plot.SRA} \linkS4class{SRA}
#' @importFrom dplyr %>%
#' @export
SRA_scope <- function(OM, data = list(), condition = c("catch", "effort"), selectivity = "logistic", s_selectivity = NULL, LWT = list(),
                      comp_like = c("multinomial", "lognormal"), ESS = c(30, 30),
                      max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE, ...) {

  dots <- list(...) # can be vul_par, s_vul_par, map_vul_par, map_s_vul_par, map_log_rec_dev, rescale
  if(!is.null(dots$maxF)) max_F <- dots$maxF

  comp_like <- match.arg(comp_like)
  condition <- match.arg(condition)

  dat_update <- update_SRA_data(data, OM, condition, dots)
  OM <- dat_update$OM
  data <- dat_update$data
  StockPars <- dat_update$StockPars
  FleetPars <- dat_update$FleetPars
  ObsPars <- dat_update$ObsPars

  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  nyears <- data$nyears
  nfleet <- data$nfleet
  nsurvey <- data$nsurvey

  OM@maxF <- max_F
  message("OM@maxF updated to ", max_F, ".")

  # Indices (by default selectivity of index is for total biomass)
  I_type2 <- suppressWarnings(as.numeric(data$I_type))
  I_type2[data$I_type == "B"] <- -1
  I_type2[data$I_type == "SSB"] <- -2
  I_type2[data$I_type == "est"] <- 0

  # No comp data
  if(is.null(data$CAA) && is.null(data$CAL)) {
    fix_sel <- TRUE
    message("No fishery length or age compositions were provided. Selectivity is fixed to values from OM.")
  } else {
    fix_sel <- FALSE
  }

  # Selectivity
  if(length(selectivity) == 1) selectivity <- rep(selectivity, nfleet)
  if(length(selectivity) < nfleet) stop("selectivity vector should be of length nfleet (", nfleet, ").", call. = FALSE)
  sel_test <- match(selectivity, c("logistic", "dome"))
  if(any(is.na(sel_test))) stop("selectivity vector should be either \"logistic\" or \"dome\".", call. = FALSE)
  sel <- ifelse(selectivity == "logistic", -1, 0)

  if(nsurvey > 0) {
    if(is.null(s_selectivity)) s_selectivity <- rep("logistic", nsurvey)
    if(length(s_selectivity) == 1) s_selectivity <- rep(s_selectivity, nsurvey)
    if(!any(I_type2 == 0)) {
      s_sel <- rep(1L, nsurvey)
    } else {
      if(length(s_selectivity) < nsurvey) stop("s_selectivity vector should be of length nsurvey (", nsurvey, ").", call. = FALSE)
      s_sel <- suppressWarnings(as.numeric(s_selectivity))
      s_sel[s_selectivity == "logistic"] <- -1
      s_sel[s_selectivity == "dome"] <- 0
      if(any(s_sel > 0 && I_type2 == 0)) {
        ind <- which(s_sel > 0 && I_type2 == 0)
        stop("Selectivity for survey ", paste0(ind, collapse = " "), " is estimated but s_selectivity should be either \"logistic\" or \"dome\".", call. = FALSE)
      }
    }

  } else {
    s_sel <- 1L
  }

  # Likelihood weights
  if(is.null(LWT$Chist)) {
    LWT$Chist <- rep(1, nfleet)
  } else if(length(LWT$Chist) == 1 && nfleet > 1) {
    LWT$Chist <- rep(LWT$Chist, nfleet)
  }
  if(length(LWT$Chist) != nfleet) stop("LWT$Chist should be a vector of length ", nfleet, ".")

  if(is.null(LWT$Index)) {
    LWT$Index <- rep(1, max(1, nsurvey))
  } else if(length(LWT$Index) == 1 && nsurvey > 1) {
    LWT$Index <- rep(LWT$Index, nsurvey)
  }
  if(length(LWT$Index) != max(1, nsurvey)) stop("LWT$Index should be a vector of length ", data$nsurvey, ".")

  if(is.null(LWT$CAA)) {
    LWT$CAA <- rep(1, nfleet)
  } else if(length(LWT$CAA) == 1 && nfleet > 1) {
    LWT$CAA <- rep(LWT$CAA, nfleet)
  }
  if(length(LWT$CAA) != nfleet) stop("LWT$CAA should be a vector of length ", nfleet, ".")

  if(is.null(LWT$CAL)) {
    LWT$CAL <- rep(1, nfleet)
  } else if(length(LWT$CAL) == 1 && nfleet > 1) {
    LWT$CAL <- rep(LWT$CAL, nfleet)
  }
  if(length(LWT$CAL) != nfleet) stop("LWT$CAL should be a vector of length ", nfleet, ".")

  if(is.null(LWT$ML)) {
    LWT$ML <- rep(1, nfleet)
  } else if(length(LWT$ML) == 1 && nfleet > 1) {
    LWT$ML <- rep(LWT$ML, nfleet)
  }
  if(length(LWT$ML) != nfleet) stop("LWT$ML should be a vector of length ", nfleet, ".")

  if(is.null(LWT$C_eq)) {
    LWT$C_eq <- rep(1, max(1, nfleet))
  } else if(length(LWT$C_eq) == 1 && nfleet > 1) {
    LWT$C_eq <- rep(LWT$C_eq, nfleet)
  }
  if(length(LWT$C_eq) != nfleet) stop("LWT$C_eq should be a vector of length ", nfleet, ".")

  if(is.null(LWT$s_CAA)) {
    LWT$s_CAA <- rep(1, max(1, nsurvey))
  } else if(length(LWT$s_CAA) == 1 && nsurvey > 1) {
    LWT$s_CAA <- rep(LWT$s_CAA, nsurvey)
  }
  if(length(LWT$s_CAA) != max(1, nsurvey)) stop("LWT$s_CAA should be a vector of length ", nsurvey, ".")

  if(is.null(LWT$s_CAL)) {
    LWT$s_CAL <- rep(1, max(1, nsurvey))
  } else if(length(LWT$s_CAL) == 1 && nsurvey > 1) {
    LWT$s_CAL <- rep(LWT$s_CAL, nsurvey)
  }
  if(length(LWT$s_CAL) != max(1, nsurvey)) stop("LWT$s_CAL should be a vector of length ", nsurvey, ".")

  data$LWT <- LWT

  # SR
  message(ifelse(OM@SRrel == 1, "Beverton-Holt", "Ricker"), " stock-recruitment relationship used.")

  # Fit model
  message("\nFitting model (", nsim, " simulations) ...")

  if(cores > 1 && !snowfall::sfIsRunning()) DLMtool::setup(as.integer(cores))
  if(snowfall::sfIsRunning()) {
    mod <- snowfall::sfClusterApplyLB(1:nsim, SRA_scope_est, data = data, I_type = I_type2, selectivity = sel, s_selectivity = s_sel,
                                      SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, comp_like = comp_like, ESS = ESS,
                                      max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                      FleetPars = FleetPars, dots = dots)
  } else {
    mod <- lapply(1:nsim, SRA_scope_est, data = data, I_type = I_type2, selectivity = sel, s_selectivity = s_sel,
                  SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, comp_like = comp_like, ESS = ESS,
                  max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                  FleetPars = FleetPars, dots = dots)
  }
  #assign('mod', mod, envir = globalenv())
  res <- lapply(mod, getElement, "report")
  conv <- vapply(res, getElement, logical(1), name = "conv")
  message(sum(conv), " out of ", nsim , " model fits converged (", 100*sum(conv)/nsim, "%).\n")
  if(sum(conv) < nsim) message("Non-converged iteration(s): ", paste(which(!conv), collapse = " "), "\n")
  if(sum(conv) < nsim && drop_nonconv) {
    message("Non-converged iterations will be removed.\n")
    keep <- conv
  } else {
    keep <- !logical(OM@nsim)
  }

  # Test for identical sims
  all_identical_sims_fn <- function() {
    vector_fn <- function(x) sum(mean(x) - x) == 0
    array_fn <- function(x) {
      x_mean <- apply(x, 2:length(dim(x)), mean)
      all(apply(x, 1, identical, x_mean))
    }
    run_test <- function(x) if(is.null(dim(x))) vector_fn(x) else array_fn(x)
    StockPars_subset <- StockPars[c("hs", "procsd", "ageM", "M_ageArray", "Linf", "Len_age", "Wt_age", "Mat_age")]
    if(!any(data$CAL > 0, na.rm = TRUE) || !any(data$s_CAL > 0, na.rm = TRUE) || !any(data$ML > 0, na.rm = TRUE)) {
      StockPars_subset <- c(StockPars_subset, StockPars["LenCV"])
    }
    S_test <- vapply(StockPars_subset, run_test, logical(1))
    if(data$nfleet == 1 && !any(data$CAL > 0, na.rm = TRUE) && !any(data$CAA > 0, na.rm = TRUE)) {
      FleetPars_subset <- StockPars[c("L5", "LFS", "Vmaxlen")]
      FleetPars_subset <- lapply(FleetPars_subset, function(x) x[nyears, ])
      F_test <- vapply(FleetPars_subset, run_test, logical(1))
    } else {
      F_test <- TRUE
    }
    if(data$nsurvey > 0 && !any(data$I_sd > 0, na.rm = TRUE)) {
      O_test <- run_test(ObsPars$Isd)
    } else {
      O_test <- TRUE
    }
    return(all(c(S_test, F_test, O_test)))
  }
  if(all_identical_sims_fn()) { # All identical sims detected
    mean_fit_output <- mod[[1]]
  } else if(mean_fit) { ### Fit to life history means if mean_fit = TRUE
    message("Generating additional model fit from mean values of parameters in the operating model...\n")
    mean_fit_output <- SRA_scope_est(data = data, I_type = I_type2, selectivity = sel, s_selectivity = s_sel,
                                     SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, comp_like = comp_like, ESS = ESS,
                                     max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars,
                                     FleetPars = FleetPars, mean_fit = TRUE, dots = dots)
  } else mean_fit_output <- list()
  if(length(mean_fit_output) > 0 && !mean_fit_output$report$conv) warning("Mean fit model did not appear to converge.")

  ### R0
  OM@cpars$R0 <- vapply(1:length(mod), function(x) ifelse("log_R0" %in% names(mod[[x]]$obj$par), res[[x]]$R0, StockPars$R0[x]), numeric(1))
  message("Range of unfished recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))

  ### Depletion and init D
  if(any(data$C_eq > 0) || any(data$E_eq > 0)) {
    initD <- vapply(res, function(x) x$E[1]/x$E0[1], numeric(1))
    message("Estimated range in initial spawning depletion: ", paste(round(range(initD), 2), collapse = " - "))
  }

  OM@cpars$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0_SR, numeric(1))
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")

  ### Selectivity and F
  ### Find
  OM@isRel <- FALSE
  if(nfleet > 1) {
    F_matrix <- lapply(res, getElement, "F_at_age")
    apical_F <- lapply(F_matrix, function(x) apply(x, 1, max))
    Find <- do.call(rbind, apical_F)

    V <- Map("/", e1 = F_matrix, e2 = apical_F)
    expand_V_matrix <- function(x) {
      y <- matrix(x[nyears, ], proyears, maxage, byrow = TRUE)
      rbind(x, y)
    }
    V <- lapply(V, expand_V_matrix)
    V2 <- array(unlist(V), c(nyears + proyears, maxage, nsim))

    OM@cpars$V <- aperm(V2, c(3, 2, 1))
    OM@cpars$Find <- Find
    message("Historical F and selectivity trends set in OM@cpars$Find and OM@cpars$V, respectively.")
    message("Selectivity during projection period is set to that in most recent historical year.")

  } else { # nfleet = 1

    OM@cpars$L5 <- vapply(res, getElement, numeric(1), "L5")
    message("Range of OM@cpars$L5: ", paste(round(range(OM@cpars$L5), 2), collapse = " - "))

    OM@cpars$LFS <- vapply(res, getElement, numeric(1), "LFS")
    message("Range of OM@cpars$LFS: ", paste(round(range(OM@cpars$LFS), 2), collapse = " - "))

    if(selectivity == "logistic") {
      OM@cpars$Vmaxlen <- rep(1, nsim)
      message("With logistic selectivity, setting OM@cpars$Vmaxlen = 1")

    } else {
      OM@cpars$Vmaxlen <- vapply(res, getElement, numeric(1), "Vmaxlen")
      message("Range of OM@cpars$Vmaxlen: ", paste(round(range(OM@cpars$Vmaxlen), 2), collapse = " - "))
    }

    OM@cpars$Find <- t(do.call(cbind, lapply(res, getElement, "F")))
    message("Historical F set in OM@cpars$Find.")
  }

  Eff <- apply(OM@cpars$Find, 2, range)
  OM@EffLower <- Eff[1, ]
  OM@EffUpper <- Eff[2, ]
  if(length(OM@EffYears) != nyears) OM@EffYears <- 1:nyears
  message("Historical effort trends set in OM@EffLower and OM@EffUpper.\n")

  ### Rec devs
  OM@cpars$Perr <- StockPars$procsd
  make_Perr <- function(x) {
    bias_corr <- ifelse(x$obj$env$data$est_rec_dev, exp(-0.5 * x$report$tau^2), 1)
    res <- exp(x$report$log_rec_dev) * bias_corr
    res[1] <- res[1] * x$report$R_eq/x$report$R0
    return(res)
  }
  Perr <- do.call(rbind, lapply(mod, make_Perr))

  make_early_Perr <- function(x) {
    res <- x$report$R_eq * x$report$NPR_equilibrium / (x$report$R0 * x$report$NPR_unfished[[1]])
    return(rev(res[-1]))
  }
  early_Perr <- do.call(rbind, lapply(mod, make_early_Perr))

  OM@cpars$Perr_y <- StockPars$Perr_y
  OM@cpars$Perr_y[, 1:(OM@maxage - 1)] <- early_Perr
  OM@cpars$Perr_y[, OM@maxage:(OM@maxage + nyears - 1)] <- Perr

  log_rec_dev <- do.call(rbind, lapply(res, getElement, "log_rec_dev"))
  OM@cpars$AC <- apply(log_rec_dev, 1, function(x) acf(x, lag.max = 1, plot = FALSE)$acf[2])
  OM@AC <- range(OM@cpars$AC)

  pro_Perr_y <- matrix(rnorm(proyears * nsim, rep(StockPars$procmu, proyears), rep(StockPars$procsd, proyears)),
                       nsim, proyears)
  for(y in 2:proyears) pro_Perr_y[, y] <- OM@cpars$AC * pro_Perr_y[, y - 1] + pro_Perr_y[, y] * sqrt(1 - OM@cpars$AC^2)
  OM@cpars$Perr_y[, (OM@maxage+nyears):ncol(OM@cpars$Perr_y)] <- exp(pro_Perr_y)

  message("Recruitment standard deviations set in OM@cpars$Perr.")
  message("Historical recruitment trends set in OM@cpars$Perr_y.")
  message("Range of recruitment autocorrelation OM@AC: ", paste(round(range(OM@AC), 2), collapse = " - "))
  message("Future recruitment deviations sampled with autocorrelation (in OM@cpars$Perr_y).\n")

  ### Assign OM variables that were used in the SRA to the output
  OM@cpars$Len_age <- StockPars$Len_age
  OM@cpars$LenCV <- StockPars$LenCV
  OM@cpars$Wt_age <- StockPars$Wt_age

  if(any(apply(StockPars$Mat_age, 1, function(x) all(x >= 0.5)))) { # Any simulations where all mat_age > 0.5?
    OM@cpars$L50 <- StockPars$L50
    OM@cpars$L95 <- StockPars$L95
  } else {
    OM@cpars$Mat_age <- StockPars$Mat_age
  }
  OM@cpars$M_ageArray <- StockPars$M_ageArray

  OM@cpars$h <- StockPars$hs

  OM@cpars$plusgroup <- rep(1L, nsim)
  OM@cpars$Iobs <- ObsPars$Iobs
  message("Growth, maturity, natural mortality, and steepness values from SRA are set in OM@cpars.\n")

  ### Output list
  E <- do.call(rbind, lapply(res[keep], getElement, "E"))
  N <- array(sapply(res[keep], getElement, "N"), c(nyears+1, maxage, sum(keep)))
  CAA_pred <- array(sapply(res[keep], getElement, "CAApred"), c(nyears, maxage, nfleet, sum(keep)))
  CAL_pred <- array(sapply(res[keep], getElement, "CALpred"), c(nyears, length(data$length_bin), nfleet, sum(keep)))

  output <- new("SRA", OM = Sub_cpars(OM, keep), SSB = E, NAA = aperm(N, c(3, 1, 2)), CAA = aperm(CAA_pred, c(4, 1:3)),
                CAL = aperm(CAL_pred, c(4, 1:3)), mean_fit = mean_fit_output, conv = conv[keep], data = data, Misc = res[keep])

  # Data in cpars
  if(sum(output@data$Chist > 0, na.rm = TRUE) || nsurvey > 0) {

    real_Data <- new("Data")
    real_Data@Year <- (output@OM@CurrentYr - output@OM@nyears + 1):output@OM@CurrentYr
    if(sum(output@data$Chist > 0, na.rm = TRUE) && all(!is.na(output@data$Chist))) {
      real_Data@Cat <- matrix(rowSums(output@data$Chist, na.rm = TRUE), 1, nyears)
      real_Data@CV_Cat <- matrix(sqrt(exp(0.01^1 - 1)), 1, nyears)
      message("Historical catch data added to OM@cpars$Data@Cat with default catch CV = 0.01.")
    }
    if(.hasSlot(real_Data, "AddInd") && nsurvey > 0) {
      real_Data@AddInd <- array(output@data$Index, c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))
      real_Data@CV_AddInd <- array(sqrt(exp(output@data$I_sd^2) - 1), c(nyears, nsurvey, output@OM@nsim)) %>%
        aperm(perm = c(3, 2, 1))

      AddIndV <- lapply(output@Misc, function(x) x$s_vul[nyears, , , drop = FALSE]) %>% unlist() %>% array(dim = c(maxage, nsurvey, OM@nsim))
      real_Data@AddIndV <- aperm(AddIndV, c(3, 2, 1))
      message("Historical indices added to OM@cpars$Data@AddInd.")
    }
    output@OM@cpars$Data <- real_Data
  }

  message("Complete.")

  return(output)
}

#### Parameters
# R0
# sel
# rec_devs

#### Data list
# Cat vector by year
# Index matrix by year and fleet
# Index_type vector by fleet
# CAA matrix by year, age
# CAL matrix by year, length bin
# mean length vector by year

#### Life history
# M array by nsim, age, year
# mat_age array by nsim, age
# Length at age array by nsim, age
# Weight at age arraty by nsim, length
# CV_LAA vector by nsim
# h steepness

#### Other settings
# SR type
# selectivity type
# maxage integer
SRA_scope_est <- function(x = 1, data, I_type, selectivity, s_selectivity, SR_type = c("BH", "Ricker"), LWT = list(),
                          comp_like = c("multinomial", "lognormal"), ESS = c(30, 30),
                          max_F = 3, integrate = FALSE, StockPars, ObsPars, FleetPars, mean_fit = FALSE,
                          control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3), dots = list()) {

  SR_type <- match.arg(SR_type)
  comp_like <- match.arg(comp_like)

  nyears <- data$nyears
  nfleet <- data$nfleet
  max_age <- dim(data$CAA)[2]
  nsurvey <- ncol(data$Index)

  data$CAA <- apply(data$CAA, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  data$CAL <- apply(data$CAL, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  CAA_n <- apply(data$CAA, c(1, 3), sum, na.rm = TRUE)
  CAL_n <- apply(data$CAL, c(1, 3), sum, na.rm = TRUE)

  data$s_CAA <- apply(data$s_CAA, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  data$s_CAL <- apply(data$s_CAL, c(1, 3), SRA_tiny_comp) %>% aperm(c(2, 1, 3))
  s_CAA_n <- apply(data$s_CAA, c(1, 3), sum, na.rm = TRUE)
  s_CAL_n <- apply(data$s_CAL, c(1, 3), sum, na.rm = TRUE)

  if(comp_like == "multinomial") {
    for(ff in 1:nfleet) { # Annual sums to effective sample size
      data$CAA[,,ff] <- data$CAA[,,ff]/CAA_n[,ff] * pmin(CAA_n[,ff], ESS[1])
      data$CAL[,,ff] <- data$CAL[,,ff]/CAL_n[,ff] * pmin(CAL_n[,ff], ESS[2])
    }
    CAA_n <- pmin(CAA_n, ESS[1])
    CAL_n <- pmin(CAL_n, ESS[2])

    for(sur in 1:nsurvey) { # Annual sums to effective sample size
      data$s_CAA[,,sur] <- data$s_CAA[,,sur]/s_CAA_n[,sur] * pmin(s_CAA_n[,sur], ESS[1])
      data$s_CAL[,,sur] <- data$s_CAL[,,sur]/s_CAL_n[,sur] * pmin(s_CAL_n[,sur], ESS[2])
    }
    s_CAA_n <- pmin(s_CAA_n, ESS[1])
    s_CAL_n <- pmin(s_CAL_n, ESS[2])
  }

  LWT_C <- matrix(c(LWT$Chist, LWT$CAA, LWT$CAL, LWT$ML, LWT$C_eq), nrow = nfleet, ncol = 5)
  LWT_Index <- cbind(LWT$Index, LWT$s_CAA, LWT$s_CAL)

  if(mean_fit) {
    # Average across simulations for arrays: M_ageArray, Len_age, Mat_age (mean across index 1)
    # Matrices: L5, LFS, Vmaxlen (mean across rows)
    # Vectors: Isd, Len_CV, hs, R0 (effort only)
    mean_vector <- function(x) rep(mean(x), length(x))
    mean_matrix <- function(x) matrix(apply(x, 1, mean), nrow(x), ncol(x))
    mean_array <- function(x) {
      means <- apply(x, c(2, 3), mean)
      new_output <- array(means, dim(x)[c(2,3,1)])
      return(aperm(new_output, c(3,1,2)))
    }

    StockPars_ind <- match(c("M_ageArray", "Len_age", "Mat_age"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_array)

    StockPars_ind <- match(c("hs", "LenCV", "procsd"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_vector)

    StockPars_ind <- match("ageM", names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_matrix)

    if(data$condition == "effort") StockPars$R0 <- mean_vector(StockPars$R0)

    FleetPars_ind <- match(c("L5", "LFS", "Vmaxlen"), names(FleetPars))
    FleetPars[FleetPars_ind] <- lapply(FleetPars[FleetPars_ind], mean_matrix)

    ObsPars$Isd <- mean_vector(ObsPars$Isd)
  }

  if(is.null(data$I_sd)) data$I_sd <- matrix(sdconv(1, ObsPars$Isd[x]), nyears, nsurvey)

  if(!is.null(data$Chist) && any(data$Chist > 0, na.rm = TRUE)) {
    if(is.null(dots$rescale)) {
      rescale <- 1/mean(data$Chist, na.rm = TRUE)
    } else {
      rescale <- dots$rescale
    }
    C_hist <- data$Chist * rescale
  } else {
    rescale <- 1
    C_hist <- matrix(0, nyears, nfleet)
  }

  if(!is.null(data$Ehist) && any(data$Ehist > 0, na.rm = TRUE)) {
    rescale_effort <- 1/mean(data$Ehist, na.rm = TRUE)
    E_hist <- data$Ehist * rescale_effort
  } else {
    rescale_effort <- 1
    E_hist <- matrix(0, nyears, nfleet)
  }

  TMB_data_all <- list(condition = data$condition,
                       nll_C = as.integer((all(CAA_n <= 0) & all(CAL_n <= 0) & all(is.na(data$ML)) & all(is.na(data$Index))) || nfleet > 1), # if condition = "effort"
                       I_hist = data$Index, sigma_I = data$I_sd, CAA_hist = data$CAA, CAA_n = pmin(CAA_n, ESS[1]),
                       CAL_hist = data$CAL, CAL_n = pmin(CAL_n, ESS[2]), s_CAA_hist = data$s_CAA, s_CAA_n = s_CAA_n,
                       s_CAL_hist = data$s_CAL, s_CAL_n = s_CAL_n, length_bin = data$length_bin, mlen = data$ML,
                       n_y = nyears, max_age = max_age, nfleet = nfleet, nsurvey = nsurvey,
                       M = t(StockPars$M_ageArray[x, , 1:nyears]), len_age = t(StockPars$Len_age[x, , 1:(nyears+1)]),
                       Linf = StockPars$Linf[x], CV_LAA = StockPars$LenCV[x], wt = t(StockPars$Wt_age[x, , 1:(nyears+1)]),
                       mat = t(StockPars$Mat_age[x, , 1:(nyears+1)]), vul_type = selectivity, s_vul_type = s_selectivity, I_type = I_type, abs_I = data$abs_I,
                       I_basis = data$I_basis, SR_type = SR_type, LWT_C = LWT_C, LWT_Index = LWT_Index, comp_like = comp_like,
                       max_F = max_F, rescale = rescale, ageM = min(nyears, ceiling(StockPars$ageM[x, 1])),
                       est_early_rec_dev = rep(0, max_age-1), yindF = as.integer(rep(0.5 * nyears, nfleet)))

  if(data$condition == "catch") {
    TMB_data <- list(model = "SRA_scope", C_hist = C_hist, C_eq = data$C_eq * rescale, E_hist = E_hist, E_eq = rep(0, nfleet))
  } else {
    TMB_data <- list(model = "SRA_scope", C_hist = C_hist, C_eq = rep(0, nfleet), E_hist = E_hist, E_eq = data$E_eq * rescale_effort)
  }

  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)

  LFS <- FleetPars$LFS[nyears, x]
  L5 <- FleetPars$L5[nyears, x]
  Vmaxlen <- FleetPars$Vmaxlen[nyears, x]

  if(is.null(dots$vul_par)) {
    vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, nfleet)
  } else {
    vul_par <- dots$vul_par
  }
  vul_par[2, ] <- log(vul_par[1, ] - vul_par[2, ])
  vul_par[1, ] <- logit(pmin(vul_par[1, ]/StockPars$Linf[x]/0.95, 0.95))
  vul_par[3, ] <- logit(pmin(vul_par[3, ], 0.99))

  if(is.null(dots$map_vul_par)) {
    map_vul_par <- matrix(0, 3, nfleet)
    map_vul_par[3, as.logical(selectivity)] <- NA

    for(ff in 1:nfleet) {
      if(all(data$CAA[,,ff] <= 0, na.rm = TRUE) && all(data$CAL[,,ff] <= 0, na.rm = TRUE)) map_vul_par[, ff] <- NA
    }
    if(any(!is.na(map_vul_par))) map_vul_par[!is.na(map_vul_par)] <- 1:sum(!is.na(map_vul_par))
  } else {
    map_vul_par <- dots$map_vul_par
  }

  # s_vul_par, and map
  if(is.null(dots$s_vul_par)) {
    s_vul_par <- matrix(c(LFS, L5, Vmaxlen), 3, nsurvey)
  } else {
    s_vul_par <- dots$s_vul_par
  }
  s_vul_par[2, ] <- log(s_vul_par[1, ] - s_vul_par[2, ])
  s_vul_par[1, ] <- logit(min(s_vul_par[1, ]/StockPars$Linf[x]/0.95, 0.95))
  s_vul_par[3, ] <- logit(pmin(s_vul_par[3, ], 0.99))

  if(is.null(dots$map_s_vul_par)) {
    map_s_vul_par <- matrix(0, 3, nsurvey)
    map_s_vul_par[3, s_selectivity < 0] <- NA # if logistic
    for(sur in 1:nsurvey) {
      if(I_type[sur] != 0 || (all(data$s_CAA[,,sur] <= 0, na.rm = TRUE) & all(data$s_CAL[,,sur] <= 0, na.rm = TRUE))) {
        map_s_vul_par[, sur] <- NA
      }
    }
    if(any(!is.na(map_s_vul_par))) map_s_vul_par[!is.na(map_s_vul_par)] <- 1:sum(!is.na(map_s_vul_par))
  } else {
    map_s_vul_par <- dots$map_s_vul_par
  }

  log_F_start <- matrix(0, nyears, nfleet)
  log_F_start[TMB_data$yindF - 1, 1:nfleet] <- log(0.1)
  TMB_params <- list(log_R0 = ifelse(TMB_data_all$nll_C, 0, log(StockPars$R0[x])),
                     transformed_h = transformed_h, vul_par = vul_par, s_vul_par = s_vul_par,
                     log_q_effort = rep(log(0.1), nfleet),
                     log_F = log_F_start, log_F_equilibrium = rep(log(0.05), nfleet),
                     log_sigma_mlen = log(data$ML_sd), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = rep(0, max_age - 1), log_rec_dev = rep(0, nyears))

  map <- list()
  if(data$condition == "effort" && !TMB_data_all$nll_C) map$log_R0 <- factor(NA)
  map$transformed_h <- map$log_tau <- factor(NA)
  map$vul_par <- factor(map_vul_par)
  map$s_vul_par <- factor(map_s_vul_par)
  if(data$condition == "catch") {
    map$log_q_effort <- factor(rep(NA, nfleet))
    if(any(data$C_eq == 0)) {
      map_log_F_equilibrium <- rep(NA, nfleet)
      map_log_F_equilibrium[data$C_eq > 0] <- 1:sum(data$C_eq > 0)
      map$log_F_equilibrium <- factor(map_log_F_equilibrium)
    }
  } else {
    map$log_F <- factor(matrix(NA, nyears, nfleet))
    map$log_F_equilibrium <- factor(rep(NA, nfleet))
  }

  map$log_sigma_mlen <- factor(rep(NA, nfleet))
  map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
  if(is.null(dots$map_log_rec_dev)) {
    map$log_rec_dev <- factor(1:nyears)
  } else {
    map$log_rec_dev <- factor(dots$map_log_rec_dev)
  }
  TMB_data$est_rec_dev <- ifelse(is.na(map$log_rec_dev), 0, 1)

  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev") else random <- NULL

  obj <- MakeADFun(data = c(TMB_data, TMB_data_all), parameters = TMB_params, map = map, random = random,
                   inner.control = inner.control, DLL = "MSEtool", silent = TRUE)

  mod <- optimize_TMB_model(obj, control, restart = 0)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)

  if(rescale != 1) {

    vars_div <- c("B", "E", "Cat", "C_eq_pred", "CAApred", "CALpred", "s_CAApred", "s_CALpred", "CN", "Cpred", "N", "N_full", "VB",
                  "R", "R_early", "R_eq", "VB0", "R0", "B0", "E0", "N0", "E0_SR")
    vars_mult <- c("Brec", "q")
    var_trans <- c("R0", "q")
    fun_trans <- c("/", "*")

    if(data$condition == "catch" || TMB_data_all$nll_C) {
      fun_fixed <- c("log", NA)
      rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
    } else if(any(data$Chist > 0, na.rm = TRUE)) {
      rescale <- 1/exp(mean(log(data$Chist/report$Cpred), na.rm = TRUE))

      fun_fixed <- c(NA, NA)
      rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
    }

  }

  report$F_at_age <- report$Z - report$M[1:nyears, ]
  report$vul_len <- get_vul_len(report)
  report <- get_s_vul_len(report, TMB_data_all)
  report$rescale <- rescale

  return(list(obj = obj, opt = opt, SD = SD, report = c(report, list(conv = !is.character(opt) && SD$pdHess))))
}


#' @rdname SRA_scope
#' @export
plot_SRA_scope <- function(...) {
  .Deprecated("report_SRA_scope", msg = "plot_SRA_scope is now deprecated in favor of plot() which generates a markdown report.")
}

#' @rdname SRA_scope
#' @export
Sub_cpars <- function(OM, sims = 1:OM@nsim) {

  if(is.numeric(sims)) {
    sims2 <- logical(OM@nsim)
    sims2[sims] <- TRUE
  } else if(is.logical(sims) && length(sims) == OM@nsim) {
    sims2 <- sims
  } else stop("Logical vector sims need to be of length ", OM@nsim)

  if(any(!sims2)) {
    message("Removing simulations: ", paste0(which(!sims2), collapse = " "))
    cpars <- OM@cpars

    subset_Data <- function(xx, Data, sims) {
      z <- slot(Data, xx)
      if(!all(is.na(z))) z <- z[sims, , , drop = FALSE]
      return(z)
    }
    subset_function <- function(x, sims) {
      if(is.matrix(x)) {
        return(x[sims, , drop = FALSE])
      } else if(is.array(x)) {
        if(length(dim(x)) == 3) return(x[sims, , , drop = FALSE])
        if(length(dim(x)) == 4) return(x[sims, , , , drop = FALSE])
        if(length(dim(x)) == 5) return(x[sims, , , , , drop = FALSE])
      } else if(class(x)[[1]] == "Data") {
        s_names <- c("AddIndV", "AddInd", "CV_AddInd")
        update_slots <- lapply(s_names, subset_Data, Data = x, sims = sims)
        for(i in 1:length(s_names)) slot(x, s_names[i]) <- update_slots[[i]]
        return(x)
      } else return(x[sims])
    }

    cpars2 <- lapply(cpars, subset_function, sims = sims2)
    OM@cpars <- cpars2
    OM@nsim <- sum(sims2)

    message("Set OM@nsim = ", OM@nsim)
  }

  return(OM)
}


get_vul_len <- function(report) {
  sls <- (report$LFS - report$L5)/sqrt(-log(0.05, 2))
  srs <- (report$Linf - report$LFS)/sqrt(-log(report$Vmaxlen, 2))

  asc <- Map(function(x, y) 2^-((report$length_bin - y)/x * (report$length_bin - y)/x), x = sls, y = report$LFS)
  dsc <- Map(function(x, y, z) ifelse(z > rep(0.99, length(report$length_bin)), 1, 2^-((report$length_bin - y)/x * (report$length_bin - y)/x)),
             x = srs, y = report$LFS, z = report$Vmaxlen)
  vul <- Map(function(x, y, z) ifelse(report$length_bin > x, y, z), x = report$LFS, y = dsc, z = asc)
  do.call(cbind, vul)
}

get_s_vul_len <- function(report, TMB_data) {
  s_vul_len <- matrix(NA, length(TMB_data$length_bin), TMB_data$nsurvey) # length-based: matrix of dimension nlbin, nsurvey

  for(i in 1:TMB_data$nsurvey) {
    if(TMB_data$I_type[i] == 0) {
      sls <- (report$s_LFS[i] - report$s_L5[i])/sqrt(-log(0.05, 2))
      srs <- (report$Linf - report$s_LFS[i])/sqrt(-log(report$s_Vmaxlen[i], 2))

      asc <- 2^-((report$length_bin - report$s_LFS[i])/sls * (report$length_bin - report$s_LFS[i])/sls)
      dsc <- ifelse(report$s_Vmaxlen[i] > rep(0.99, length(report$length_bin)), 1,
                    2^-((report$length_bin - report$s_LFS[i])/srs * (report$length_bin - report$s_LFS[i])/srs))
      s_vul_len[, i] <- ifelse(report$length_bin > report$s_LFS[i], dsc, asc)
    }
    if(TMB_data$I_type[i] > 0) s_vul_len[, i] <- report$vul_len[, TMB_data$I_type[i]]
  }

  report$s_vul_len <- s_vul_len
  return(report)
}
