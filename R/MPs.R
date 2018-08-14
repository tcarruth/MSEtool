#
#' Data-rich management procedures
#'
#' A suite of data-rich management procedures (MPs) included in the package. Additional MPs,
#' with specific model configurations (e.g., stock-recruit function or fixing certain parameters) or alternative
#' ramped harvest control rules can be created with \link{make_MP} and the available Assess and HCR objects.
#'
#' @name Data-rich-MP
#' @param x A position in the Data object.
#' @param Data An object of class Data
#' @param reps Numeric, the number of stochastic replicates for the management advice.
#'
#' @examples
#' \dontrun{
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "SCA_MSY", "SCA_4010"))
#' }
#' @return An object of class Rec which contains the management recommendation.
NULL





#' @describeIn Data-rich-MP A statistical catch-at-age model with a TAC recommendation based on fishing at UMSY,
#' and default arguments for configuring \link{SCA}.
#' @export
SCA_MSY <- function(x, Data, reps = 1, ...) {
  do_Assessment <- SCA(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps)
  return(Rec)
}
class(SCA_MSY) <- "MP"
#SCA_MSY <- make_MP(SCA, HCR_MSY)


#' @describeIn Data-rich-MP An SCA with a TAC recommendation based on fishing at 75\% of UMSY.
#' @export
SCA_75MSY <- function(x, Data, reps = 1, ...) {
  do_Assessment <- SCA(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps, MSY_frac = 0.75)
  return(Rec)
}
class(SCA_75MSY) <- "MP"
#SCA_75MSY <- make_MP(SCA, HCR_MSY, MSY_frac = 0.75)


#' @describeIn Data-rich-MP An SCA with a 40-10 control rule.
#' @export
SCA_4010 <- function(x, Data, reps = 1, ...) {
  do_Assessment <- SCA(x = x, Data = Data, ...)
  Rec <- HCR40_10(do_Assessment, reps = reps)
  return(Rec)
}
class(SCA_4010) <- "MP"
#SCA_4010 <- make_MP(SCA, HCR40_10)



#' @describeIn Data-rich-MP A state-space delay difference model with a TAC recommendation based on fishing at UMSY,
#' and default arguments for configuring \link{DD_SS}.
#' @export
DDSS_MSY <- function(x, Data, reps = 1, ...) {
  do_Assessment <- DD_SS(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps)
  return(Rec)
}
class(DDSS_MSY) <- "MP"
#DDSS_MSY <- make_MP(DD_SS, HCR_MSY)


#' @describeIn Data-rich-MP A state-space delay difference model with a TAC recommendation based on fishing at 75\% of UMSY.
#' @export
DDSS_75MSY <- function(x, Data, reps = 1, ...) {
  do_Assessment <- DD_SS(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps, MSY_frac = 0.75)
  return(Rec)
}
class(DDSS_75MSY) <- "MP"
#DDSS_75MSY <- make_MP(DD_SS, HCR_MSY, MSY_frac = 0.75)


#' @describeIn Data-rich-MP A state-space delay difference model with a 40-10 control rule.
#' @export
DDSS_4010 <- function(x, Data, reps = 1, ...) {
  do_Assessment <- DD_SS(x = x, Data = Data, ...)
  Rec <- HCR40_10(do_Assessment, reps = reps)
  return(Rec)
}
class(DDSS_4010) <- "MP"
#DDSS_4010 <- make_MP(DD_SS, HCR40_10)



#' @describeIn Data-rich-MP A surplus production model with a TAC recommendation based on fishing at UMSY,
#' and default arguments for configuring \link{SP}.
#' @export
SP_MSY <- function(x, Data, reps = 1, ...) {
  do_Assessment <- SP(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps)
  return(Rec)
}
class(SP_MSY) <- "MP"
#SP_MSY <- make_MP(SP, HCR_MSY)


#' @describeIn Data-rich-MP A surplus production model with a TAC recommendation based on fishing at 75\% of UMSY.
#' @export
SP_75MSY <- function(x, Data, reps = 1, ...) {
  do_Assessment <- SP(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps, MSY_frac = 0.75)
  return(Rec)
}
class(SP_75MSY) <- "MP"
#SP_75MSY <- make_MP(SP, HCR_MSY, MSY_frac = 0.75)


#' @describeIn Data-rich-MP A surplus production model with a 40-10 control rule.
#' @export
SP_4010 <- function(x, Data, reps = 1, ...) {
  do_Assessment <- SP(x = x, Data = Data, ...)
  Rec <- HCR40_10(do_Assessment, reps = reps)
  return(Rec)
}
class(SP_4010) <- "MP"
#SP_4010 <- make_MP(SP, HCR40_10)
