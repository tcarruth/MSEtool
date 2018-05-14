#' Make a custom management procedure (MP)
#'
#' Function operator that combines a function of class \code{Assess} and a function of
#' class \code{HCR} to create a management procedure (MP). The resulting function
#' can then be tested in closed-loop simulation via \code{\link[DLMtool]{runMSE}}.
#'
#' @param .Assess A function of class \code{Assess}.
#' @param .HCR A function of class \code{HCR}.
#' @param diagnostic (TRUE/FALSE) Whether assessment output will be written to
#' \code{\link[DLMtool]{DLMenv}}. Intended for use with \code{\link[DLMtool]{runMSE}}.
#' @param ... Additional arguments to be passed to \code{.Assess} and \code{.HCR}.
#'
#' @return A function of class \code{MP}.
#'
#' @examples
#' # A delay-difference model with a 40-10 control rule
#' DD_40_10 <- make_MP(DD_TMB, HCR40_10)
#'
#' \dontrun{
#' myOM <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "DD_40_10"))
#'
#' ls(DLMenv) # Model output during MSE is assigned to this environment.
#' plot_diagnostic()
#' }
#'
#' # MP that uses a surplus production model which assumes Binit_frac = 0.5.
#' SP_MSY <- make_MP(SP, HCR_MSY, Binit_frac = 0.5)
#' @importFrom pryr make_function
#' @seealso \link{plot_diagnostics}
#' @export
make_MP <- function(.Assess, .HCR, diagnostic = FALSE, ...) {
  if(is.character(.Assess)) {
    .Assess <- as.symbol(.Assess)
  } else {
    .Assess <- substitute(.Assess)
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

  if(diagnostic) {
    MP_body <- bquote({
      opt_timing <- system.time(do_Assessment <- .(.Assess)(x = x, Data = Data, ...))[3]
      Rec <- .(.HCR)(do_Assessment, reps = reps, ...)
      Assess_diagnostic()
      return(Rec)
    })
  } else {
    MP_body <- bquote({
      do_Assessment <- .(.Assess)(x = x, Data = Data, ...)
      Rec <- .(.HCR)(do_Assessment, reps = reps, ...)
      return(Rec)
    })
  }

  custom_MP <- make_function(args = alist(x = , Data = , reps = 1),
                             body = MP_body)
  formals(custom_MP) <- c(formals(custom_MP), list(...))
  class(custom_MP) <- "MP"
  return(custom_MP)
}


