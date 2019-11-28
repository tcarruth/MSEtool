

SRA_scope_data <- function (OM, Data, ...) {

  # < code for making Chist, Index, ML, CAA, CAL, I_sd matrices out of list of Data objects >
  nyears<-OM@nyears
  CYcond<-length(Data@Cat[1,]) != nyears
  if(CYcond) message(paste0("Catch data for a different duration than operating model. Data@Cat[1,] is of length ",length(Data@Cat[1,]),", but OM@nyears is ",nyears))
  if(CYcond)stop("OM, Data, incompatibility")

  Chist<-Data@Cat[1, ]
  Index<-Data@Ind[1, ]
  I_sd<-Data@CV_Ind[1, ]
  CAA<-Data@CAA[1, , ]
  CAA<-NULL
  CAL<-Data@CAL[1, , ]
  CAL<-NULL
  ML<-Data@ML[1, ]
  ML<-NULL
  length_bin<-Data@CAL_mids

  out<-SRA_scope(OM, Chist=Chist, Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL,
                            ML = ML, length_bin = length_bin, ...)
  return(out)
}



update_SRA_data <- function(data, OM, condition, dots) {

  assign_for_compatibility <- function(x) {
    if(is.null(data[[x]]) && !is.null(dots[[x]])) {
      data[[x]] <<- getElement(dots, x)
      dots[[x]] <<- NULL
    }
    invisible()
  }

  dat_names <- c("Chist", "Ehist", "Index", "I_sd", "I_type", "CAA", "CAL", "ML", "ML_sd", "C_eq", "E_eq", "s_CAA", "s_CAL", "length_bin")

  lapply(dat_names, assign_for_compatibility)


  if(is.null(data$Chist) && !is.null(data$Ehist)) condition <- "effort"

  if(condition == "catch") {
    if(is.null(data$Chist)) {
      stop("Full time series of catch is needed.", call. = FALSE)
    } else {
      if(any(is.na(data$Chist))) {
        stop("One or more of the historical annual catch observations is missing. Suggestion: use linear interpolation to fill these data.", call. = FALSE)
      }
      if(any(data$Chist < 0)) stop("All catch values should be zero or greater.", call. = FALSE)

      # Convert single fleet inputs to multiple fleet, e.g. matrices to arrays
      if(!is.matrix(data$Chist)) data$Chist <- matrix(data$Chist, ncol = 1)

      data$nyears <- nrow(data$Chist)
      data$nfleet <- ncol(data$Chist)
      data$Ehist <- matrix(0, data$nyears, data$nfleet)

      if(is.null(data$Index) && is.null(data$CAA) && is.null(data$CAL) && is.null(data$ML) && is.null(data$Ehist)) {
        message("No data other than Chist is provided. Model will switch to conditioning on equilibrium effort.")
        data$condition <- "effort"
        data$Ehist <- matrix(1, data$nyears, data$nfleet)
      } else {
        data$condition <- "catch"
      }
    }
  }

  if(condition == "effort") {
    data$condition <- "effort"
    if(is.null(data$Ehist)) {
      stop("Full time series of effort is needed.")
    } else {
      if(any(is.na(data$Ehist))) stop("Effort time series is not complete (contains NA's")
      if(any(data$Ehist < 0)) stop("All effort values should be positive.")

      if(!is.matrix(data$Ehist)) data$Ehist <- matrix(data$Ehist, ncol = 1)

      data$nyears <- nrow(data$Ehist)
      data$nfleet <- ncol(data$Ehist)

      if(!is.null(data$Chist) && !is.matrix(data$Chist)) data$Chist <- matrix(data$Chist, ncol = 1)
      if(is.null(data$Chist)) data$Chist <- matrix(0, data$nyears, data$nfleet)
    }
  }

  message("SRA model is conditioned on ", data$condition)
  message(data$nfleet, " fleet(s) detected.")
  message(data$nyears, " years of data detected.")

  # Match number of historical years of catch/effort to OM
  if(OM@nyears != data$nyears) {
    cpars_cond <- length(OM@cpars) > 0 && any(vapply(OM@cpars, function(x) class(x) == "matrix" || class(x) == "array", logical(1)))
    if(cpars_cond) {
      stmt <- paste0("OM@nyears != length(", ifelse(data$condition == "catch", "Chist", "Ehist"), "). ",
                     "There will be indexing errors in your custom parameters (OM@cpars).")
      stop(stmt, call. = FALSE)
    } else {
      message("OM@nyears was updated to length(", ifelse(data$condition == "catch", "Chist", "Ehist"), "): ", data$nyears)
      OM@nyears <- data$nyears
    }
  }
  if(length(OM@CurrentYr) == 0) OM@CurrentYr <- data$nyears

  # Indices
  if(!is.null(data$Index)) {
    if(is.vector(data$Index)) {
      if(length(data$Index) != data$nyears) stop("Length of Index vector does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
      data$Index <- matrix(data$Index, ncol = 1)
    } else if(is.matrix(data$Index)) {
      if(nrow(data$Index) != data$nyears) stop("Number of rows of Index matrix does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    } else stop("Index is neither a vector nor a matrix.", call. = FALSE)

    data$nsurvey <- ncol(data$Index)

    # Match index to I_type
    if(is.null(data$I_type)) data$I_type <- rep("B", data$nsurvey)
    if(length(data$I_type) != ncol(data$Index)) {
      stop("Length of I_type needs to be ", data$nsurvey, call. = FALSE)
    }

    I_type_check <- match(data$I_type, c("B", "SSB", "est", 1:data$nfleet))
    if(any(is.na(I_type_check))) stop("I_type vector needs to be entries of either: \"est\", \"SSB\", \"B\", or 1 - ", data$nfleet, ".", call. = FALSE)
  } else {
    data$nsurvey <- 0
    data$Index <- matrix(NA, ncol = 1, nrow = data$nyears)
    data$I_type <- "B"
  }

  if(!is.null(data$I_sd)) {
    if(is.vector(data$I_sd)) {
      if(length(data$I_sd) != data$nyears) stop("Length of I_sd vector does not equal nyears (", data$nyears, ").", call. = FALSE)
      data$I_sd <- matrix(data$I_sd, ncol = 1)
    } else if(is.matrix(data$I_sd)) {
      if(nrow(data$I_sd) != data$nyears) stop("Number of rows of I_sd matrix does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
      if(ncol(data$I_sd) != data$nsurvey) stop("Number of columns of I_sd matrix does not equal nsurvey (", data$nsurvey, ").", call. = FALSE)
    }
  }

  message(data$nsurvey, " survey(s) detected.")

  # Process age comps
  if(!is.null(data$CAA)) {

    if(is.matrix(data$CAA)) data$CAA <- array(data$CAA, c(dim(data$CAA), 1))

    if(dim(data$CAA)[1] != data$nyears) {
      stop("Number of CAA rows (", dim(data$CAA)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$CAA)[2] < OM@maxage) {
      message("Number of CAA columns (", dim(data$CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no observations for ages greater than ", dim(data$CAA)[2], " and filling with zeros.")
      add_ages <- OM@maxage - dim(data$CAA)[2]
      CAA_new <- array(0, c(data$nyears, OM@maxage, data$nfleet))
      CAA_new[, 1:dim(data$CAA)[2], ] <- data$CAA
      data$CAA <- CAA_new
    }
    if(dim(data$CAA)[2] > OM@maxage) {
      OM@maxage <- dim(data$CAA)[2]
      message("Increasing OM@maxage to ", OM@maxage, ".")
    }
    if(dim(data$CAA)[3] != data$nfleet) {
      stop("Number of CAA slices (", dim(data$CAA)[3], ") does not equal nfleet (", data$nfleet, "). NAs are acceptable.", call. = FALSE)
    }

  } else {
    data$CAA <- array(0, c(data$nyears, OM@maxage, data$nfleet))
  }

  # Sample life history, selectivity, and obs parameters
  set.seed(OM@seed)
  StockPars <- SampleStockPars(OM, msg = FALSE)
  ObsPars <- SampleObsPars(OM)
  FleetPars <- SampleFleetPars(OM, msg = FALSE)

  # Process length comps
  if(!is.null(data$CAL)) {
    if(is.matrix(data$CAL)) data$CAL <- array(data$CAL, c(dim(data$CAL), 1))

    if(is.null(data$length_bin)) {
      stop("You must specify length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(data$CAL)[1] != data$nyears) {
      stop("Number of CAL rows (", dim(data$CAL)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$CAL)[2] != length(data$length_bin)) {
      stop("Number of CAL columns (", dim(data$CAL)[2], ") does not equal length(length_bin) (", length(data$length_bin), ").", call. = FALSE)
    }
    if(dim(data$CAL)[3] != data$nfleet) {
      stop("Number of CAL slices (", dim(data$CAA)[3], ") does not equal nfleet (", data$nfleet, "). NAs are acceptable.", call. = FALSE)
    }
  } else {
    data$CAL <- array(0, c(data$nyears, length(StockPars$CAL_binsmid), data$nfleet))
    data$length_bin <- StockPars$CAL_binsmid
  }

  # Process mean lengths
  if(!is.null(data$ML)) {
    if(is.vector(data$ML)) {
      if(length(data$ML) != data$nyears) stop("Mean length vector (ML) must be of length ", data$nyears, ".", call. = FALSE)
      data$ML <- matrix(data$ML, ncol = 1)
    }
    if(nrow(data$ML) != data$nyears) stop("Number of ML rows (", nrow(data$ML), ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    if(ncol(data$ML) != data$nfleet) stop("Number of ML columns (", ncol(data$ML), ") does not equal nfleet (", data$nfleet, "). NAs are acceptable.", call. = FALSE)

    if(is.null(data$ML_sd)) {
      data$ML_sd <- apply(data$ML, 2, mean, na.rm = TRUE)
    } else if(length(data$ML_sd) == 1) data$ML_sd <- rep(data$ML_sd, data$nfleet)
    if(length(data$ML_sd) != data$nfleet) stop("Mean length SD vector (ML_sd) must be of length ", data$nfleet, ".", call. = FALSE)
  } else {
    data$ML <- matrix(NA, nrow = data$nyears, ncol = data$nfleet)
    data$ML_sd <- rep(0.1, data$nfleet)
  }

  # Process equilibrium catch/effort - Ceq
  if(is.null(data$C_eq)) data$C_eq <- rep(0, data$nfleet)
  if(data$condition == "catch") {
    if(length(data$C_eq) == 1) data$C_eq <- rep(data$C_eq, data$nfleet)
    if(length(data$C_eq) < data$nfleet) stop("C_eq needs to be of length nfleet (", data$nfleet, ").", call. = FALSE)
  }
  if(is.null(data$E_eq)) data$E_eq <- rep(0, data$nfleet)
  if(data$condition == "effort") {
    if(length(data$E_eq) == 1) data$E_eq <- rep(data$E_eq, data$nfleet)
    if(length(data$E_eq) < data$nfleet) stop("E_eq needs to be of length nfleet (", data$nfleet, ").", call. = FALSE)
  }

  # Process survey age comps
  if(!is.null(data$s_CAA)) {

    if(is.matrix(data$s_CAA)) data$s_CAA <- array(data$s_CAA, c(dim(data$s_CAA), 1))

    if(dim(data$s_CAA)[1] != data$nyears) {
      stop("Number of s_CAA rows (", dim(data$s_CAA)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$s_CAA)[2] < OM@maxage) {
      message("Number of s_CAA columns (", dim(data$s_CAA)[2], ") does not equal OM@maxage (", OM@maxage, ").")
      message("Assuming no observations for ages greater than ", dim(data$s_CAA)[2], " and filling with zeros.")
      add_ages <- OM@maxage - dim(data$s_CAA)[2]
      CAA_new <- array(0, c(data$nyears, OM@maxage, data$survey))
      CAA_new[, 1:dim(data$s_CAA)[2], ] <- data$s_CAA
      data$s_CAA <- CAA_new
    }
    if(dim(data$s_CAA)[2] > OM@maxage) {
      stop("Error in age dimension of s_CAA.", call. = FALSE)
    }
    if(dim(data$s_CAA)[3] != data$nsurvey) {
      stop("Number of CAA slices (", dim(data$s_CAA)[3], ") does not equal nsurvey (", data$nsurvey, "). NAs are acceptable.", call. = FALSE)
    }

  } else {
    data$s_CAA <- array(0, c(data$nyears, OM@maxage, ncol(data$Index)))
  }

  # Process survey length comps
  if(!is.null(data$s_CAL)) {
    if(is.matrix(data$s_CAL)) data$s_CAL <- array(data$s_CAL, c(dim(data$s_CAL), 1))

    if(is.null(data$length_bin)) {
      stop("You must specify length_bin, which is the mean length of each length bin (columns) of the CAL data.", call. = FALSE)
    }
    if(dim(data$s_CAL)[1] != data$nyears) {
      stop("Number of s_CAL rows (", dim(data$s_CAL)[1], ") does not equal nyears (", data$nyears, "). NAs are acceptable.", call. = FALSE)
    }
    if(dim(data$s_CAL)[2] != length(data$length_bin)) {
      stop("Number of s_CAL columns (", dim(data$s_CAL)[2], ") does not equal length(length_bin) (", length(data$length_bin), ").", call. = FALSE)
    }
    if(dim(data$s_CAL)[3] != data$nsurvey) {
      stop("Number of s_CAL slices (", dim(data$s_CAL)[3], ") does not equal nsurvey (", data$nsurvey, "). NAs are acceptable.", call. = FALSE)
    }
  } else {
    data$s_CAL <- array(0, c(data$nyears, length(data$length_bin), ncol(data$Index)))
  }

  return(list(data = data, OM = OM, StockPars = StockPars, ObsPars = ObsPars, FleetPars = FleetPars))

}
