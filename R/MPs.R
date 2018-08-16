
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
#' @param ... Additional arguments passed to the Assessment model.
#'
#' @examples
#' \dontrun{
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "SCA_MSY", "SCA_4010"))
#' }
#' @return An object of class Rec which contains the management recommendation.
NULL




#temp <- lapply(avail("Assess"), function(x) paste(format(match.fun(x)), collapse = " "))
#Assess_dep <- lapply(temp, DLMtool:::match_slots)
#names(Assess_dep) <- avail("Assess")
Assess_dep <- list(DD_SS = "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge",
                   DD_TMB = "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge",
                   SCA = "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge",
                   SCA2 = "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge",
                   SP = "Data@Cat, Data@Ind",
                   SP_SS = "Data@Cat, Data@Ind, Data@CV_Ind")


get_dependencies <- function(Assess, arg = list()) {
  dep <- getElement(Assess_dep, Assess)
  formals2 <- formals(Assess)
  arg_match <- pmatch(names(arg), names(formals2))
  formals2[arg_match[!is.na(arg_match)]] <- arg[!is.na(arg_match)]

  dep_match <- pmatch(names(dep_args), names(formals2)) # From formals
  more_deps <- dep_args[do.call(c, formals2[dep_match])]
  if(Assess == "DD_SS" && any(names(more_deps) == "fix_sigma")) more_deps$fix_sigma <- "Data@CV_Cat"
  more_deps <- paste(do.call(c, more_deps), collapse = ", ")
  paste(c(dep, more_deps), collapse = ", ")
}

dep_args <- list(fix_h = "Data@steep", fix_sigma = "Data@CV_Ind", fix_tau = "Data@sigmaR")




#' @describeIn Data-rich-MP A statistical catch-at-age model with a TAC recommendation based on fishing at UMSY,
#' and default arguments for configuring \link{SCA}.
#' @export
SCA_MSY <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$SCA)
  do_Assessment <- SCA(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps)
  return(Rec)
}))
class(SCA_MSY) <- "MP"
#SCA_MSY <- make_MP(SCA, HCR_MSY)


#' @describeIn Data-rich-MP An SCA with a TAC recommendation based on fishing at 75\% of UMSY.
#' @export
SCA_75MSY <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$SCA)
  do_Assessment <- SCA(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps, MSY_frac = 0.75)
  return(Rec)
}))
class(SCA_75MSY) <- "MP"
#SCA_75MSY <- make_MP(SCA, HCR_MSY, MSY_frac = 0.75)


#' @describeIn Data-rich-MP An SCA with a 40-10 control rule.
#' @export
SCA_4010 <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$SCA)
  do_Assessment <- SCA(x = x, Data = Data, ...)
  Rec <- HCR40_10(do_Assessment, reps = reps)
  return(Rec)
}))
class(SCA_4010) <- "MP"
#SCA_4010 <- make_MP(SCA, HCR40_10)



#' @describeIn Data-rich-MP A state-space delay difference model with a TAC recommendation based on fishing at UMSY,
#' and default arguments for configuring \link{DD_SS}.
#' @export
DDSS_MSY <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$DD_SS)
  do_Assessment <- DD_SS(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps)
  return(Rec)
}))
class(DDSS_MSY) <- "MP"
#DDSS_MSY <- make_MP(DD_SS, HCR_MSY)


#' @describeIn Data-rich-MP A state-space delay difference model with a TAC recommendation based on fishing at 75\% of UMSY.
#' @export
DDSS_75MSY <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$DD_SS)
  do_Assessment <- DD_SS(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps, MSY_frac = 0.75)
  return(Rec)
}))
class(DDSS_75MSY) <- "MP"
#DDSS_75MSY <- make_MP(DD_SS, HCR_MSY, MSY_frac = 0.75)


#' @describeIn Data-rich-MP A state-space delay difference model with a 40-10 control rule.
#' @export
DDSS_4010 <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$DD_SS)
  do_Assessment <- DD_SS(x = x, Data = Data, ...)
  Rec <- HCR40_10(do_Assessment, reps = reps)
  return(Rec)
}))
class(DDSS_4010) <- "MP"
#DDSS_4010 <- make_MP(DD_SS, HCR40_10)



#' @describeIn Data-rich-MP A surplus production model with a TAC recommendation based on fishing at UMSY,
#' and default arguments for configuring \link{SP}.
#' @export
SP_MSY <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$SP)
  do_Assessment <- SP(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps)
  return(Rec)
}))
class(SP_MSY) <- "MP"
#SP_MSY <- make_MP(SP, HCR_MSY)


#' @describeIn Data-rich-MP A surplus production model with a TAC recommendation based on fishing at 75\% of UMSY.
#' @export
SP_75MSY <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$SP)
  do_Assessment <- SP(x = x, Data = Data, ...)
  Rec <- HCR_MSY(do_Assessment, reps = reps, MSY_frac = 0.75)
  return(Rec)
}))
class(SP_75MSY) <- "MP"
#SP_75MSY <- make_MP(SP, HCR_MSY, MSY_frac = 0.75)


#' @describeIn Data-rich-MP A surplus production model with a 40-10 control rule.
#' @export
SP_4010 <- eval(bquote(function(x, Data, reps = 1, ...) {
  dependencies <- .(Assess_dep$SP)
  do_Assessment <- SP(x = x, Data = Data, ...)
  Rec <- HCR40_10(do_Assessment, reps = reps)
  return(Rec)
}))
class(SP_4010) <- "MP"
#SP_4010 <- make_MP(SP, HCR40_10)






#' Make a custom management procedure (MP)
#'
#' Function operator that combines a function of class \code{Assess} and a function of
#' class \code{HCR} to create a management procedure (MP). The resulting function
#' can then be tested in closed-loop simulation via \code{\link[DLMtool]{runMSE}}.
#'
#' @param .Assess A function of class \code{Assess}.
#' @param .HCR A function of class \code{HCR}.
#' @param diagnostic A character string describing if any additional diagnostic information from the
#' assessment models will be collected during a call with \code{\link[DLMtool]{runMSE}} (\code{"none"} is the
#' default). \code{"min"} (minimal) will collect information on convergence and \code{"full"} will also collect the
#' Assessment object generated by the \code{.Assess}. This information will be written to \code{\link[DLMtool]{DLMenv}}.
#' See example.
#' @param ... Additional arguments to be passed to \code{.Assess} and \code{.HCR}.
#'
#' @return A function of class \code{MP}.
#'
#' @examples
#' # A delay-difference model with a 40-10 control rule
#' DD_40_10 <- make_MP(DD_TMB, HCR40_10)
#'
#' # A delay difference model that will produce convergence diagnostics
#' DD_40_10 <- make_MP(DD_TMB, HCR40_10, diagnostic = "min")
#'
#' # MP that uses a Delay-Difference which assumes a Ricker stock-recruit function.
#' DD_Ricker <- make_MP(DD_TMB, HCR_MSY, SR = "Ricker")
#'
#' \dontrun{
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "DD_40_10"))
#'
#' ls(DLMenv) # Model output during MSE is assigned to this environment.
#' diagnostic_AM(myMSE)
#' }
#'
#' @importFrom pryr make_function
#' @seealso \link{diagnostic_AM} \link{retrospective_AM}
#' @export
make_MP <- function(.Assess, .HCR, diagnostic = c("none", "min", "full"), ...) {
  diagnostic <- match.arg(diagnostic)
  if(is.character(.Assess)) {
    Assess_char <- .Assess
    .Assess <- as.symbol(.Assess)
  } else {
    .Assess <- substitute(.Assess)
    Assess_char <- as.character(.Assess)
  }
  if(is.character(.HCR)) {
    .HCR <- as.symbol(.HCR)
  } else {
    .HCR <- substitute(.HCR)
  }
  if(!inherits(eval(.Assess), "Assess")) {
    stop(paste(.Assess, "does not belong to class 'Assess'. Use: avail('Assess') to find eligible objects."))
  }
  if(!inherits(eval(.HCR), "HCR")) {
    stop(paste(.HCR, "does not belong to class 'HCR.' Use: avail('HCR') to find eligible objects."))
  }

  if(diagnostic == "none") {
    MP_body <- bquote({
      dependencies <- .(get_dependencies(Assess_char, list(...)))
      do_Assessment <- .(.Assess)(x = x, Data = Data, ...)
      Rec <- .(.HCR)(do_Assessment, reps = reps, ...)
      return(Rec)
    })
  } else {
    MP_body <- bquote({
      dependencies <- .(get_dependencies(Assess_char, list(...)))
      do_Assessment <- .(.Assess)(x = x, Data = Data, ...)
      Rec <- .(.HCR)(do_Assessment, reps = reps, ...)
      Assess_diagnostic(include_assessment = .(diagnostic == "full"))
      return(Rec)
    })
  }

  custom_MP <- make_function(args = alist(x = , Data = , reps = 1), body = MP_body)
  formals(custom_MP) <- c(formals(custom_MP), list(...))
  class(custom_MP) <- "MP"
  return(custom_MP)
}
