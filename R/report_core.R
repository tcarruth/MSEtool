

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
profile_likelihood <- function(Assessment, figure = TRUE, save_figure = TRUE,
                               save_dir = tempdir(), ...) {
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
retrospective <- function(Assessment, nyr = 5, figure = TRUE, save_figure = TRUE,
                          save_dir = tempdir()) {
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


