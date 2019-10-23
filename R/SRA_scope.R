


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
#' @param Chist A vector of historical catch, should be of length OM@@nyears. If there are multiple fleets: a matrix of OM@@nyears rows and nfleet columns.
#' Ideally, the first year of the catch series represents unfished conditions (see also \code{C_eq}).
#' @param Ehist A vector of historical effort, should be of length OM@@nyears (see also \code{E_eq}).
#' @param condition String to indicate whether the SRA model is conditioned on catch or effort.
#' @param Index A vector of values of an index (of length OM@@nyears). If there are multiple surveys: a matrix of historical indices of abundances, with rows
#' indexing years and columns indexing surveys.
#' @param I_sd A vector or matrix of standard deviations (lognormal distribution) for the indices corresponding to the entries in \code{Index}.
#' If NULL, this function will use values from OM@@Iobs.
#' @param CAA Age composition matrix with nyears rows and OM@@maxage columns. If multiple fleets: an array with dimension: nyears, OM@@maxage, and nfleets.
#' @param CAL Length composition matrix with nyears rows and columns indexing the length bin. If multiple fleets: an array with dimension: nyears,
#' length bins, and nfleets.
#' @param ML A vector of mean length observations (length OM@@nyears), or if multiple fleets: matrix of dimension: nyears and nfleets. Generally, should not
#' be used if \code{CAL} is also provided, unless mean length and length comps are independently sampled.
#' @param length_bin A vector for the midpoints of the length bins for \code{CAL}. All length bin widths should be equal in size.
#' @param C_eq A numeric vector of length nfleet for the equilibrium catch for each fleet in \code{Chist} prior to the first year of the operating model.
#' Zero implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.
#' @param E_eq The equilibrium effort for each fleet in \code{Ehist} prior to the first year of the operating model.
#' Zero implies unfished conditions in year one. Otherwise, this is used to estimate depletion in the first year of the data.
#' @param ML_sd The standard deviation (normal distribution) of the observed mean lengths. If there are multiple fleets, a vector of length nfleet.
#' If \code{NULL}, default value is \code{0.1 * mean(ML)}.
#' @param selectivity A character vector of length nfleet to indicate \code{"logistic"} or \code{"dome"} selectivity for each fleet in \code{Chist}.
#' @param fix_dome Logical, if \code{selectivity = "dome"}, determines whether the descending limb of selectivity is fixed or not.
#' @param I_type A character vector (length nsurvey) to indicate the type of biomass for which each index follows. Either \code{"B"} for
#' total biomass, or \code{"SSB"} for spawning biomass. If \code{NULL}, "B" is used. Use numbers if the index corresponds to a fleet in \code{Chist}.
#' @param LWT A named list of likelihood weights for the SRA model. See details.
#' @param ESS A numeric vector of length two for the maximum effective samples size of the age and length compositions, respectively for the
#' multinomial likelihood function. The annual sample size of an age or length composition sample is the minimum of ESS or the number of observations.
#' @param max_F The maximum F for any fleet in the scoping model (higher F's in the model are penalized in the objective function).
#' @param cores Integer for the number of CPU cores for the stock reduction analysis.
#' @param integrate Logical, whether to treat recruitment deviations as penalized parameters (FALSE) or random effects (TRUE).
#' @param mean_fit Logical, whether to run an additional with mean values of life history parameters from the OM.
#' @param sims A logical vector of length \code{OM@@nsim} or a numberic vector indicating which simulations to keep.
#' @param drop_nonconv Logical, whether to drop non-converged fits of the SRA model.
#' @param ... Other arguments to add in the future.
#' @return An object of class \linkS4class{SRA} (see link for output dimensions).
#'
#' @details
#' One of indices, age compositions, or length compositions should be provided in addition to the historical catch.
#' Selectivity is fixed to values sampled from \code{OM} if no age or length compositions are provided.
#'
#' \code{LWT} is a named list containing the likelihood weights (values > 0) with the possible options:
#' \itemize{
#' \item Chist: A vector of length nfleet.
#' \item Index: A vector of length nsurvey.
#' \item CAA, CAL, ML, C_eq: A vector of length nfleet for each.
#' }
#' By default, all likelihood weights are equal to one if not specified by the user. Likelihoods for CAA and CAL can also be adjusted by changing the
#' multinomial sample size. See argument \code{ESS}.
#'
#' \code{plot_SRA_scope} is now deprecated in favor of \link{plot.SRA}
#'
#' Output from \code{SRA_scope} is placed in objects in \code{OM@@cpars}. \code{Sub_cpars} is a convenient function to subset simulations
#' for the operating model, for example, to remove simulations from unconverged model fits or outlier simulations.
#'
#' @note If the operating model \code{OM} uses time-varying growth or M, then those trends will be used in the SRA as well.
#' Time-varying life history parameters create ambiguity in the calculation and interpretation of depletion and reference points in \link[DLMtool]{runMSE}.
#' See section D.5 of \code{DLMtool::userguide()}.
#'
#' Here, the initial depletion (OM@@cpars$initD) is calculated based on unfished spawning biomass using growth and M in the first year. If growth/M varies,
#' then this reference point may no longer be relevant at the end of the historical period and in the projection period.
#'
#' The easiest way to turn off time-varying growth/M is by setting: \code{OM@@Msd <- OM@@Linfsd <- OM@@Ksd <- c(0, 0)}.
#'
#' @author Q. Huynh
#' @seealso \link{plot.SRA} \linkS4class{SRA}
#' @export
SRA_scope <- function(OM, Chist = NULL, Ehist = NULL, condition = c("catch", "effort"), Index = NULL, I_sd = NULL, CAA = NULL, CAL = NULL, ML = NULL, length_bin = NULL,
                      C_eq = 0, E_eq = 0, ML_sd = NULL, selectivity = "logistic", I_type = NULL, LWT = list(), ESS = c(30, 30),
                      fix_dome = FALSE, max_F = 3, cores = 1L, integrate = FALSE, mean_fit = FALSE, drop_nonconv = FALSE, ...) {

  condition <- match.arg(condition)
  if(is.null(Chist) && !is.null(Ehist)) condition <- "effort"
  message("SRA model is conditioned on ", condition)

  if(condition == "catch") {
    if(is.null(Chist)) {
      stop("Full time series of catch is needed.")
    } else {
      if(any(is.na(Chist))) {
        stop("One or more of the historical annual catch observations is missing. Suggestion: use linear interpolation to fill these data.")
      }
      if(any(Chist < 0)) stop("All catch values should be zero or greater.")

      # Convert single fleet inputs to multiple fleet, e.g. matrices to arrays
      if(!is.matrix(Chist)) Chist <- matrix(Chist, ncol = 1)

      nyears <- nrow(Chist)
      nfleet <- ncol(Chist)
      Ehist <- matrix(0, nyears, nfleet)

      if(is.null(Index) && is.null(CAA) && is.null(CAL) && is.null(ML) && is.null(Ehist)) {
        condition <- "effort"
        message("No data other than Chist is provided. Model will switch to conditioning on equilibrium effort.")
        Ehist <- matrix(1, nyears, nfleet)
      }
    }
  }

  if(condition == "effort") {
    if(is.null(Ehist)) {
      stop("Full time series of effort is needed.")
    } else {
      if(any(is.na(Ehist))) stop("Effort time series is not complete (contains NA's")
      if(any(Ehist < 0)) stop("All effort values should be positive.")

      if(!is.matrix(Ehist)) Ehist <- matrix(Ehist, ncol = 1)

      nyears <- nrow(Ehist)
      nfleet <- ncol(Ehist)

      if(!is.null(Chist) && !is.matrix(Chist)) Chist <- matrix(Chist, ncol = 1)
      if(is.null(Chist)) Chist <- matrix(0, nyears, nfleet)
    }
  }

  nsim <- OM@nsim
  proyears <- OM@proyears
  maxage <- OM@maxage
  message(nfleet, " fleet(s) detected.")

  OM@maxF <- max_F
  message("OM@maxF updated to ", max_F, ".")

  # Match number of historical years of catch to OM
  if(OM@nyears != nyears) {
    cpars_cond <- length(OM@cpars) > 0 && any(vapply(OM@cpars, function(x) class(x) == "matrix" || class(x) == "array", logical(1)))
    stmt <- ifelse(cpars_cond, ". This could create indexing errors in your custom parameters (OM@cpars).", "")
    warning("OM@nyears will be updated to length(Chist): ", nyears, stmt)
    OM@nyears <- nyears
  }
  if(length(OM@CurrentYr) == 0) OM@CurrentYr <- nyears

  # Indices
  if(!is.null(Index)) {
    if(is.vector(Index)) {
      if(length(Index) != nyears) stop("Length of Index vector does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
      Index <- matrix(Index, ncol = 1)
    } else if(is.matrix(Index)) {
      if(nrow(Index) != nyears) stop("Number of rows of Index matrix does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    } else stop("Index is neither a vector nor a matrix.", call. = FALSE)

    nsurvey <- ncol(Index)

    # Match index to I_type
    if(is.null(I_type)) I_type <- rep("B", ncol(Index))
    if(length(I_type) != ncol(Index)) {
      stop("Length of I_type needs to be ", ncol(Index), call. = FALSE)
    }

    I_type_check <- match(I_type, c("B", "SSB", 1:nfleet))
    if(any(is.na(I_type_check))) stop("I_type vector needs to be entries of either: \"SSB\", \"B\", or 1 - ", ncol(Chist), ".", call. = FALSE)
  } else {
    nsurvey <- 0
    Index <- matrix(NA, ncol = 1, nrow = nyears)
    I_type <- "B"
  }
  message(nsurvey, " survey(s) detected.")

  I_type2 <- suppressWarnings(as.numeric(I_type))
  I_type2[I_type == "B"] <- -1
  I_type2[I_type == "SSB"] <- -2

  if(!is.null(I_sd)) {
    if(is.vector(I_sd)) {
      if(length(I_sd) != nyears) stop("Length of I_sd vector does not equal nyears (", nyears, ").", call. = FALSE)
      I_sd <- matrix(I_sd, ncol = 1)
    } else if(is.matrix(I_sd)) {
      if(nrow(I_sd) != nyears) stop("Number of rows of I_sd matrix does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(I_sd) != nsurvey) stop("Number of columns of I_sd matrix does not equal nsurvey (", nsurvey, ").", call. = FALSE)
    } else I_sd <- NULL
  }

  # No comp data
  if(is.null(CAA) && is.null(CAL)) {
    fix_sel <- TRUE
    message("No length or age compositions were provided. Selectivity is fixed to values from OM.")
  } else {
    fix_sel <- FALSE
  }

  # Selectivity
  if(!fix_sel) {
    sel_test <- match(selectivity, c("logistic", "dome"))
    if(any(is.na(sel_test))) stop("selectivity vector should be either \"logistic\" or \"dome\".", call. = FALSE)

    if(length(selectivity) < nfleet) stop("selectivity vector should be of length nfleet (", nfleet, ").", call. = FALSE)
  }
  sel <- ifelse(selectivity == "logistic", 1L, 0L)

  # Process age comps
  if(!is.null(CAA)) {

    if(is.matrix(CAA)) CAA <- array(CAA, c(dim(CAA), 1))

    if(dim(CAA)[1] != nyears) {
      stop("Number of CAA rows (", dim(CAA)[1], ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(CAA)[2] < OM@maxage) {
      message("Number of CAA columns (", dim(CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no CAA for ages greater than ", dim(CAA)[2], " and filling with zeros.")
      add_ages <- OM@maxage - dim(CAA)[2]
      CAA_new <- array(0, c(nyears, OM@maxage, nfleet))
      CAA_new[, 1:dim(CAA)[2], ] <- CAA
      CAA <- CAA_new
    }
    if(dim(CAA)[2] > OM@maxage) {
      OM@maxage <- maxage <- dim(CAA)[2]
      message("Increasing OM@maxage to ", OM@maxage, ".")
    }
    if(dim(CAA)[3] != nfleet) {
      stop("Number of CAA slices (", dim(CAA)[3], ") does not equal nfleet (", nfleet, "). NAs are acceptable.", call. = FALSE)
    }

  } else {
    CAA <- array(0, c(nyears, maxage, nfleet))
  }

  # Sample life history, selectivity, and obs parameters
  set.seed(OM@seed)
  StockPars <- SampleStockPars(OM, msg = FALSE)
  ObsPars <- SampleObsPars(OM)
  FleetPars <- SampleFleetPars(OM, msg = FALSE)

  # Process length comps
  if(!is.null(CAL)) {
    if(is.matrix(CAL)) {
      CAL2 <- array(NA, c(dim(CAL), 1))
      CAL2[,,1] <- CAL
      CAL <- CAL2
    }
    if(is.null(length_bin)) {
      stop("You must specify the argument length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(CAL)[1] != nyears) {
      stop("Number of CAL rows (", dim(CAL)[1], ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(CAL)[2] != length(length_bin)) {
      stop("Number of CAL columns (", dim(CAL)[2], ") does not equal length(length_bin) (", length(length_bin), ").", call. = FALSE)
    }
    if(dim(CAL)[3] != nfleet) {
      stop("Number of CAL slices (", dim(CAA)[3], ") does not equal nfleet (", nfleet, "). NAs are acceptable.", call. = FALSE)
    }
  } else {
    CAL <- array(0, c(nyears, length(StockPars$CAL_binsmid), nfleet))
    length_bin <- StockPars$CAL_binsmid
  }

  # Process mean lengths
  if(!is.null(ML)) {
    if(is.vector(ML)) {
      if(length(ML) != nyears) stop("Mean length vector (ML) must be of length ", nyears, ".", call. = FALSE)
      ML <- matrix(ML, ncol = 1)
    }
    if(nrow(ML) != nyears) stop("Number of ML rows (", nrow(ML), ") does not equal nyears (", nyears, "). NAs are acceptable.", call. = FALSE)
    if(ncol(ML) != nfleet) stop("Number of ML columns (", ncol(ML), ") does not equal nfleet (", nfleet, "). NAs are acceptable.", call. = FALSE)

    if(is.null(ML_sd)) {
      ML_sd <- apply(ML, 2, mean, na.rm = TRUE)
    } else if(length(ML_sd) == 1) ML_sd <- rep(ML_sd, nfleet)
    if(length(ML_sd) != nfleet) stop("Mean length SD vector (ML_sd) must be of length ", nfleet, ".", call. = FALSE)

  } else {
    ML <- matrix(NA, nrow = nyears, ncol = nfleet)
    ML_sd <- rep(0.1, nfleet)
  }

  # Process equilibrium catch - Ceq
  if(condition == "catch") {
    if(length(C_eq) == 1) C_eq <- rep(C_eq, nfleet)
    if(length(C_eq) < nfleet) stop("C_eq needs to be of length nfleet (", nfleet, ").", call. = FALSE)
  }
  if(condition == "effort") {
    if(length(E_eq) == 1) E_eq <- rep(E_eq, nfleet)
    if(length(E_eq) < nfleet) stop("E_eq needs to be of length nfleet (", nfleet, ").", call. = FALSE)
  }

  #### Data summary
  SRA_data <- list(Chist = Chist, Ehist = Ehist, Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, ML = ML, length_bin = length_bin,
                   C_eq = C_eq, E_eq = E_eq, nfleet = nfleet, nsurvey = nsurvey)

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
  if(length(LWT$Index) != max(1, nsurvey)) stop("LWT$Index should be a vector of length ", nsurvey, ".")

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

  # SR
  message(ifelse(OM@SRrel == 1, "Beverton-Holt", "Ricker"), " stock-recruitment relationship used.")

  # Fit model
  message("\nFitting model (", nsim, " simulations) ...")

  if(cores > 1 && !snowfall::sfIsRunning()) DLMtool::setup(as.integer(cores))
  if(snowfall::sfIsRunning()) {
    mod <- snowfall::sfClusterApplyLB(1:nsim, SRA_scope_est, Catch = Chist, Effort = Ehist, condition = condition,
                              Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, ML = ML,
                              ML_sd = ML_sd, length_bin = length_bin,
                              I_type = I_type2, C_eq = C_eq, E_eq = E_eq, selectivity = sel,
                              fix_selectivity = fix_sel, fix_dome = fix_dome, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, ESS = ESS,
                              max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)

  } else {
    mod <- lapply(1:nsim, SRA_scope_est, Catch = Chist, Effort = Ehist, condition = condition,
                  Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, ML = ML,
                  ML_sd = ML_sd, length_bin = length_bin,
                  I_type = I_type2, C_eq = C_eq, E_eq = E_eq, selectivity = sel,
                  fix_selectivity = fix_sel, fix_dome = fix_dome, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"), LWT = LWT, ESS = ESS,
                  max_F = max_F, integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)
  }
  #assign('mod', mod, envir = globalenv())
  res <- lapply(mod, getElement, "report")
  conv <- vapply(res, getElement, logical(1), name = "conv")
  message(sum(conv), " out of ", nsim , " model fits converged (", 100*sum(conv)/nsim, "%).\n")
  if(sum(conv) < nsim) message("Non-converged iteration(s): ", paste(which(!conv), collapse = " "), "\n")
  if(sum(conv) < nsim && drop_nonconv) {
    message("Non-converged iterations will be removed. Setting OM@nsim = ", sum(conv), "\n")
    OM@nsim <- sum(conv)
    keep <- conv
  } else {
    keep <- !logical(OM@nsim)
  }

  ### Fit to life history means if mean_fit = TRUE
  if(mean_fit) {
    message("Generating additional model fit from mean values of parameters in the operating model...\n")
    mean_fit_output <- SRA_scope_est(Catch = Chist, Effort = Ehist, condition = condition,
                                     Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, ML = ML,
                                     ML_sd = ML_sd, length_bin = length_bin,
                                     I_type = I_type2, C_eq = C_eq, E_eq = E_eq, selectivity = sel,
                                     fix_selectivity = fix_sel, fix_dome = fix_dome, SR_type = ifelse(OM@SRrel == 1, "BH", "Ricker"),
                                     LWT = LWT, ESS = ESS, max_F = max_F, mean_fit = TRUE,
                                     integrate = integrate, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars)
  }

  ### R0
  OM@cpars$R0 <- vapply(1:length(mod), function(x) ifelse("log_R0" %in% names(mod[[x]]$obj$par), res[[x]]$R0, StockPars$R0[x]), numeric(1))
  message("Range of unfished recruitment (OM@cpars$R0): ", paste(round(range(OM@cpars$R0), 2), collapse = " - "))

  ### Depletion and init D
  if(any(C_eq > 0) || any(E_eq > 0)) {
    initD <- vapply(res, function(x) x$E[1]/x$E0[1], numeric(1))
    message("Estimated range in initial spawning depletion: ", paste(round(range(initD), 2), collapse = " - "))
  }

  OM@cpars$D <- vapply(res, function(x) x$E[length(x$E)-1]/x$E0_SR, numeric(1))
  message("Range of spawning depletion (OM@cpars$D): ", paste(round(range(OM@cpars$D), 2), collapse = " - "), "\n")

  ### Selectivity and F
  ### Find
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
    V2 <- array(unlist(V), c(nyears+proyears, maxage, nsim))

    OM@cpars$V <- aperm(V2, c(3, 2, 1))
    OM@cpars$Find <- Find
    message("Historical F and selectivity trends set in OM@cpars$Find and OM@cpars$V, respectively.")
    message("Selectivity during projection period is equal to that in most recent historical year.")

  } else { # nfleet = 1

    if(!fix_sel) {
      OM@isRel <- FALSE

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
    }

    OM@cpars$Find <- t(do.call(cbind, lapply(res, getElement, "F")))
    message("Historical F set in OM@cpars$Find.")
  }

  if(fix_sel) {
    OM@cpars$L5 <- if(is.matrix(FleetPars$L5)) FleetPars$L5[nyears, ] else FleetPars$L5
    OM@cpars$LFS <- if(is.matrix(FleetPars$LFS)) FleetPars$LFS[nyears, ] else FleetPars$LFS
    OM@cpars$Vmaxlen <- if(is.matrix(FleetPars$Vmaxlen)) FleetPars$Vmaxlen[nyears, ] else FleetPars$Vmaxlen
    V <- FleetPars$V
    maxV <- apply(FleetPars$V, c(1, 3), max)
    for(i in 1:maxage) V[,i,] <- V[,i,]/maxV
    OM@cpars$V <- V
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
                       OM@nsim, proyears)
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
  CAL_pred <- array(sapply(res[keep], getElement, "CALpred"), c(nyears, length(length_bin), nfleet, sum(keep)))

  output <- new("SRA", OM = Sub_cpars(OM, keep), SSB = E, NAA = aperm(N, c(3, 1, 2)), CAA = aperm(CAA_pred, c(4, 1:3)),
                CAL = aperm(CAL_pred, c(4, 1:3)), conv = conv[keep], data = SRA_data, Misc = res[keep])
  if(mean_fit) output@mean_fit <- mean_fit_output

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

get_vul_len <- function(report) {
  sls <- (report$LFS - report$L5)/sqrt(-log(0.05, 2))
  srs <- (report$Linf - report$LFS)/sqrt(-log(report$Vmaxlen, 2))

  asc <- Map(function(x, y) 2^-((report$length_bin - y)/x * (report$length_bin - y)/x), x = sls, y = report$LFS)
  dsc <- Map(function(x, y, z) ifelse(z > rep(0.99, length(report$length_bin)), 1, 2^-((report$length_bin - y)/x * (report$length_bin - y)/x)),
             x = srs, y = report$LFS, z = report$Vmaxlen)
  vul <- Map(function(x, y, z) ifelse(report$length_bin > x, y, z), x = report$LFS, y = dsc, z = asc)
  do.call(cbind, vul)
}

SRA_scope_est <- function(x = 1, Catch = NULL, Effort = NULL, Index = NULL, condition = c("catch", "effort"),
                          I_sd = NULL, CAA = NULL, CAL = NULL, ML = NULL, ML_sd = NULL, length_bin,
                          I_type, C_eq = 0, E_eq = 0, SR_type = c("BH", "Ricker"), LWT = list(), ESS = c(30, 30),
                          StockPars, ObsPars, FleetPars, integrate = FALSE, selectivity, fix_selectivity = TRUE,
                          fix_dome = FALSE, mean_fit = FALSE, max_F = 3,
                          control = list(iter.max = 2e+05, eval.max = 4e+05), inner.control = list(maxit = 1e3)) {
  condition <- match.arg(condition)
  SR_type <- match.arg(SR_type)

  nyears <- nrow(Catch)
  if(is.null(nyears)) nyears <- nrow(Effort)
  nfleet <- ncol(Catch)
  if(is.null(nfleet)) nfleet <- ncol(Effort)
  max_age <- dim(CAA)[2]
  nsurvey <- ncol(Index)

  CAA_n <- t(rowSums(aperm(CAA, c(3, 1, 2)), dims = 2, na.rm = TRUE))
  CAL_n <- t(rowSums(aperm(CAL, c(3, 1, 2)), dims = 2, na.rm = TRUE))

  for(ff in 1:nfleet) {
    CAA[,,ff] <- CAA[,,ff]/CAA_n[,ff] * pmin(CAA_n[,ff], ESS[1])
    CAL[,,ff] <- CAL[,,ff]/CAL_n[,ff] * pmin(CAL_n[,ff], ESS[2])
  }

  LWT_C <- matrix(c(LWT$Chist, LWT$CAA, LWT$CAL, LWT$ML, LWT$C_eq), nrow = nfleet, ncol = 5)

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

    StockPars_ind <- match(c("hs", "LenCV", "procsd", "ageM"), names(StockPars))
    StockPars[StockPars_ind] <- lapply(StockPars[StockPars_ind], mean_vector)

    if(condition == "effort") StockPars$R0 <- mean_vector(StockPars$R0)

    FleetPars_ind <- match(c("L5", "LFS", "Vmaxlen"), names(FleetPars))
    FleetPars[FleetPars_ind] <- lapply(FleetPars[FleetPars_ind], mean_matrix)

    ObsPars$Isd <- mean_vector(ObsPars$Isd)
  }

  if(is.null(I_sd)) I_sd <- matrix(sdconv(1, ObsPars$Isd[x]), nyears, nsurvey)

  TMB_data_all <- list(condition = condition,
                       nll_C = as.integer((all(CAA_n <= 0) & all(CAL_n <= 0) & all(is.na(ML)) & all(is.na(Index))) || nfleet > 1), # if condition = "effort"
                       I_hist = Index, sigma_I = I_sd, CAA_hist = CAA, CAA_n = pmin(CAA_n, ESS[1]),
                       CAL_hist = CAL, CAL_n = pmin(CAL_n, ESS[2]), length_bin = length_bin, mlen = ML,
                       n_y = nyears, max_age = ncol(CAA), nfleet = nfleet, nsurvey = nsurvey,
                       M = t(StockPars$M_ageArray[x, , 1:nyears]), len_age = t(StockPars$Len_age[x, , 1:(nyears+1)]),
                       Linf = StockPars$Linf[x], CV_LAA = StockPars$LenCV[x], wt = t(StockPars$Wt_age[x, , 1:(nyears+1)]),
                       mat = t(StockPars$Mat_age[x, , 1:(nyears+1)]), vul_type = selectivity, I_type = I_type, SR_type = SR_type,
                       LWT_C = LWT_C, LWT_Index = LWT$Index, max_F = max_F, ageM = min(nyears, ceiling(StockPars$ageM[x, 1])),
                       est_early_rec_dev = rep(0, max_age-1), est_rec_dev = c(rep(1, nyears-1), 0))
  TMB_data_all$CAA_hist[TMB_data_all$CAA_hist < 1e-8] <- 1e-8
  TMB_data_all$CAL_hist[TMB_data_all$CAL_hist < 1e-8] <- 1e-8

  if(!is.null(Catch) && any(Catch > 0, na.rm = TRUE)) {
    rescale <- 1/mean(Catch, na.rm = TRUE)
    C_hist <- Catch * rescale
  } else {
    rescale <- 1
    C_hist <- matrix(0, nyears, nfleet)
  }

  if(!is.null(Effort) && any(Effort > 0, na.rm = TRUE)) {
    rescale_effort <- 1/mean(Effort, na.rm = TRUE)
    E_hist <- Effort * rescale_effort
  } else {
    rescale_effort <- 1
    E_hist <- matrix(0, nyears, nfleet)
  }

  if(condition == "catch") {
    TMB_data <- list(model = "SRA_scope", C_hist = C_hist, C_eq = C_eq * rescale,
                     E_hist = E_hist, E_eq = rep(0, nfleet))
  } else {
    TMB_data <- list(model = "SRA_scope", C_hist = C_hist, C_eq = rep(0, nfleet),
                     E_hist = E_hist, E_eq = E_eq * rescale_effort)
  }

  if(SR_type == "BH") {
    transformed_h <- logit((StockPars$hs[x] - 0.2)/0.8)
  } else transformed_h <- log(StockPars$hs[x] - 0.2)

  LFS <- FleetPars$LFS[nyears, x]
  L5 <- FleetPars$L5[nyears, x]

  vul_par <- matrix(0, 4, nfleet)
  map_vul_par <- matrix(ifelse(fix_selectivity, NA, 0), 4, nfleet)

  for(ff in 1:nfleet) {
    vul_par[1:2, ff] <- c(logit(min(LFS/StockPars$Linf[x]/0.95, 0.95)), log(LFS - L5))

    if(selectivity[ff]) { #logistic
      map_vul_par[3:4, ff] <- NA
    } else {
      vul_par[3:4, ff] <- c(-20, logit(FleetPars$Vmaxlen[nyears, x]))
      map_vul_par[3, ff] <- NA
      if(fix_dome) map_vul_par[4, ff] <- NA
    }
  }

  if(!fix_selectivity) map_vul_par[!is.na(map_vul_par)] <- 1:sum(!is.na(map_vul_par))

  TMB_params <- list(log_R0 = ifelse(TMB_data_all$nll_C, 3, log(StockPars$R0[x])),
                     transformed_h = transformed_h, vul_par = vul_par, log_q_effort = rep(log(0.1), nfleet),
                     log_F = matrix(log(0.05), nyears, nfleet), log_F_equilibrium = rep(log(0.05), nfleet),
                     log_sigma_mlen = log(ML_sd), log_tau = log(StockPars$procsd[x]),
                     log_early_rec_dev = rep(0, max_age - 1), log_rec_dev = rep(0, nyears))

  map <- list()
  if(condition == "effort" && !TMB_data_all$nll_C) map$log_R0 <- factor(NA)
  map$transformed_h <- map$log_tau <- factor(NA)
  map$vul_par <- factor(map_vul_par)
  if(condition == "catch") {
    map$log_q_effort <- factor(rep(NA, nfleet))
    if(any(C_eq == 0)) {
      map_log_F_equilibrium <- rep(NA, nfleet)
      map_log_F_equilibrium[C_eq > 0] <- 1:sum(C_eq > 0)
      map$log_F_equilibrium <- factor(map_log_F_equilibrium)
    }
  } else {
    map$log_F <- factor(matrix(NA, nyears, nfleet))
    map$log_F_equilibrium <- factor(rep(NA, nfleet))
  }

  map$log_sigma_mlen <- factor(rep(NA, nfleet))
  map$log_early_rec_dev <- factor(rep(NA, max_age - 1))
  map$log_rec_dev <- factor(c(1:(nyears-1), NA))

  if(integrate) random <- c("log_early_rec_dev", "log_rec_dev") else random <- NULL

  obj <- MakeADFun(data = c(TMB_data, TMB_data_all), parameters = TMB_params, map = map, random = random,
                   inner.control = inner.control, DLL = "MSEtool", silent = TRUE)

  mod <- optimize_TMB_model(obj, control, restart = 0)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)
  report$C_hist <- TMB_data$C_hist/rescale
  report$E_hist <- TMB_data$E_hist/rescale_effort

  vars_div <- c("B", "E", "Cat", "C_eq_pred", "CAApred", "CALpred", "CN", "Cpred", "N", "N_full", "VB",
                "R", "R_early", "R_eq", "VB0", "R0", "B0", "E0", "N0", "E0_SR")
  vars_mult <- c("Brec", "q")
  var_trans <- c("R0", "q")
  fun_trans <- c("/", "*")

  if(condition == "catch" || TMB_data_all$nll_C) {
    fun_fixed <- c("log", NA)
    rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  } else if(any(Catch > 0, na.rm = TRUE)) {
    rescale <- 1/exp(mean(log(Catch/report$Cpred), na.rm = TRUE))

    fun_fixed <- c(NA, NA)
    rescale_report(vars_div, vars_mult, var_trans, fun_trans, fun_fixed)
  }
  report$F_at_age <- report$Z - report$M[1:nyears, ]
  report$vul_len <- get_vul_len(report)
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
    subset_function <- function(x, sims) {
      if(is.matrix(x)) {
        return(x[sims, , drop = FALSE])
      } else if(is.array(x)) {
        if(length(dim(x)) == 3) return(x[sims, , , drop = FALSE])
        if(length(dim(x)) == 4) return(x[sims, , , , drop = FALSE])
        if(length(dim(x)) == 5) return(x[sims, , , , , drop = FALSE])
      } else return(x[sims])
    }

    cpars2 <- lapply(cpars, subset_function, sims = sims2)
    OM@cpars <- cpars2
    OM@nsim <- sum(sims2)

    message("Set OM@nsim = ", OM@nsim)
  }

  return(OM)
}
