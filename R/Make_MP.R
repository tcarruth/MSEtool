#' Make a custom management procedure (MP)
#'
#' Function operator that combines a function of class \code{Assess} and a function of
#' class \code{HCR} to create a management procedure (MP). The resulting function
#' can then be tested in closed-loop simulation via \code{\link[DLMtool]{runMSE}}.
#'
#' @param .Assess A function of class \code{Assess}.
#' @param .HCR A function of class \code{HCR}.
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
#' }
#'
#' # MP that uses a surplus production model which assumes Binit_frac = 0.5.
#' SP_MSY <- make_MP(SP, HCR_MSY, Binit_frac = 0.5)
#' @export make_MP
#' @importFrom pryr make_function
make_MP <- function(.Assess, .HCR, ...) {
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
  if(!inherits(eval(.Assess), "Assess")) stop(paste(.Assess, "does not belong to class 'Assess.'"))
  if(!inherits(eval(.HCR), "HCR")) stop(paste(.HCR, "does not belong to class 'HCR.'"))

  MP_body <- bquote({
    do_Assessment <- .(.Assess)(x = x, Data = Data, ...)
    Rec <- .(.HCR)(do_Assessment, reps = reps, ...)
    return(Rec)
  })

  custom_MP <- make_function(args = alist(x = , Data = , reps = 1),
                             body = MP_body)
  formals(custom_MP) <- c(formals(custom_MP), list(...))
  class(custom_MP) <- "MP"
  return(custom_MP)
}


