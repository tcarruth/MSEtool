

#' Profile likelihood of assessment models
#'
#' Profile the likelihood for leading parameters of assessment models.
#'
#' @param Assessment An S4 object of class \linkS4class{Assessment}.
#' @param figure Indicates whether a figure will be plotted.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @param ... A sequence of values of the parameter(s) for the profile.
#' See details for name of arguments to be passed on.
#' @details For the following assessment models, the required sequence of values are:
#' \itemize{
#' \item \code{DD_TMB} and \code{DD_SS}: \code{R0} and \code{h}
#' \item \code{SP} and \code{SP_SS}: \code{UMSY} and \code{MSY}
#' \item \code{SCA}: \code{R0} and \code{h}
#' \item \code{SCA2}: \code{meanR}
#' }
#' @author Q. Huynh
#' @return A data frame of negative log-likelihood values from the profile and, optionally,
#' a figure of the likelihood surface.
#' @examples
#' \donttest{
#' output <- DD_TMB(Data = DLMtool::Red_snapper)
#' pro <- profile_likelihood(output, R0 = seq(0.75, 1.25, 0.025), h = seq(0.9, 0.99, 0.01),
#' save_figure = FALSE)
#'
#' # Ensure your grid is of proper resolution. A grid that is too coarse will distort the shape of
#' # the likelihood surface.
#' }
#' @export
profile_likelihood <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  if(figure) {
    old.warning <- options()$warn
    options(warn = -1)
    on.exit(options(warn = old.warning))

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
  }

  f <- get(paste0('profile_likelihood_', Assessment@Model))
  f(Assessment, figure = figure, save_figure = save_figure, save_dir = save_dir, ...)
}


#' Retrospective analysis of assessment models
#'
#' Perform a retrospective analysis, successive removals of most recent years of data to evaluate resulting
#' parameter estimates.
#'
#' @param Assessment An S4 object of class \linkS4class{Assessment}.
#' @param nyr The maximum number of years to remove for the retrospective analysis.
#' @param figure Indicates whether plots will be drawn.
#' @return A list with an array of model output and of model estimates from
#' the retrospective analysis.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved.
#' @author Q. Huynh
#' @return Figures showing the time series of biomass and exploitation and parameter estimates
#' with successive number of years removed. For a variety of time series output (SSB, recruitment, etc.) and
#' estimates (R0, steepness, etc.), also returns a matrix of Mohn's rho (Mohn 1999).
#' @examples
#' \donttest{
#' output <- DD_TMB(Data = DLMtool::Red_snapper)
#' get_retro <- retrospective(output, nyr = 5, figure = FALSE)
#' }
#' @references
#' Mohn, R. 1999. The retrospective problem in sequential population analysis: an investigation using cod fishery
#' and simulated data. ICES Journal of Marine Science 56:473-488.
#' @export
retrospective <- function(Assessment, nyr = 5, figure = TRUE, save_figure = TRUE, save_dir = tempdir()) {
  if(figure) {
    old.warning <- options()$warn
    options(warn = -1)
    on.exit(options(warn = old.warning))

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
  }

  f <- get(paste0('retrospective_', Assessment@Model))
  f(Assessment, nyr, figure = figure, save_figure = save_figure, save_dir = save_dir)
}

#' Compare output from several assessment models
#'
#' Plot biomass, recruitment, and fishing mortality time series from several . This function can be used to compare outputs among
#' different assessment models from the same Data object.
#'
#' @param ... Objects of class \linkS4class{Assessment}.
#' @param label A character vector of the models for the legend.
#' @param color A vector of colors for each assessment model.
#' @author Q. Huynh
#' @return A set of figures of biomass, recruitment, and fishing mortality estimates among the models.
#' @examples
#' res <- cDD_SS(Data = DLMtool::SimulatedData)
#' res2 <- SCA(Data = DLMtool::SimulatedData)
#' res3 <- SCA2(Data = DLMtool::SimulatedData)
#' res4 <- VPA(Data = DLMtool::SimulatedData)
#'
#' compare_models(res, res2, res3)
#' @importFrom gplots rich.colors
#' @export
compare_models <- function(..., label = NULL, color = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  dots <- list(...)
  class_dots <- vapply(dots, inherits, logical(1), what = "Assessment")
  if(!all(class_dots)) stop("Some objects provided were not of class Assessment", call. = FALSE)

  n_assess <- length(dots)
  if(n_assess <= 1) stop("Need more than one assessment model for this function.", call. = FALSE)

  if(is.null(label) || length(label) != n_assess) label <- vapply(dots, slot, character(1), name = "Model")
  if(is.null(color) || length(color) != n_assess) {
    color <- rich.colors(n_assess)
  }

  par(mfrow = c(3, 2), mar = c(5, 4, 1, 1), oma = c(2, 0, 0, 0))

  # B/BMSY
  B_BMSY <- do.call(rbind, lapply(dots, slot, name = "SSB_SSBMSY"))
  ts_matplot(B_BMSY, expression(SSB/SSB[MSY]), color = color, dotted_one = TRUE)

  # B/B0
  B_B0 <- do.call(rbind, lapply(dots, slot, name = "SSB_SSB0"))
  ts_matplot(B_B0, expression(SSB/SSB[0]), color = color)

  # VB
  VB <- do.call(rbind, lapply(dots, slot, name = "VB"))
  ts_matplot(VB, "Vulnerable Biomass", color = color)

  # R
  RR <- lapply(dots, slot, name = "R")
  R <- match_R_years(RR)
  ts_matplot(R, "Recruitment", color = color)

  # F
  FM <- do.call(rbind, lapply(dots, slot, name = "FMort"))
  ts_matplot(FM, "Fishing Mortality", color = color)

  # F/FMSY
  F_FMSY <- do.call(rbind, lapply(dots, slot, name = "F_FMSY"))
  ts_matplot(F_FMSY, expression(F/F[MSY]), color = color, dotted_one = TRUE)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", label, col = color, xpd = TRUE, horiz = TRUE, bty = "n", lwd = 2)

  invisible()
}

#' @importFrom graphics matlines
ts_matplot <- function(m, ylab, color, dotted_one = FALSE) {
  m <- t(m)
  x <- matrix(as.numeric(rownames(m)), ncol = ncol(m), nrow = nrow(m))
  plot(NULL, NULL, xlim = range(as.numeric(rownames(m))), ylim = c(0, 1.1 * max(m, na.rm = TRUE)), xlab = "Year", ylab = ylab)
  abline(h = 0, col = "grey")
  matlines(x, m, typ = "l", col = color, lty = 1, lwd = 2)
  if(dotted_one) abline(h = 1, lty = 3)
}

match_R_years <- function(RR) {
  yrs <- do.call(c, lapply(RR, function(x) as.numeric(names(x))))
  yrs <- range(yrs)

  R <- matrix(NA, nrow = length(RR), ncol = diff(yrs) + 1)
  R_yrs <- seq(min(yrs), max(yrs))
  for(i in 1:length(RR)) {
    ind <- match(as.numeric(names(RR[[i]])), R_yrs)
    R[i, ind] <- RR[[i]]
  }
  colnames(R) <- R_yrs
  return(R)
}
