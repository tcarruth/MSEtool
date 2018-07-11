

#' Perform profile likelihood
#'
#' Performs a profile likelihood.
#'
#' @param Assessment An S4 object of class Assessment.
#' @param figure Indicates whether a figure will be plotted.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @param ... A sequence of values of the parameter(s) for the profile.
#' See details for name of arguments to be passed on.
#' @details For \code{DD_TMB}, sequence of values for \code{UMSY} and \code{MSY} are needed.
#' @author Q. Huynh
#' @return A data frame of negative log-likelihood values from profile and, optionally,
#' a figure of the likelihood surface.
#' @examples
#' \dontrun{
#' output <- DD_TMB(Data = DLMtool::Red_snapper)
#' pro <- profile_likelihood(output, R0 = seq(0.75, 1.25, 0.025), h = seq(0.9, 0.99, 0.01),
#' save_figure = FALSE)
#'
#' # Ensure your grid is of proper resolution. A grid that is too coarse will distort the shape of
#' # the likelihood surface.
#' }
#' @export
profile_likelihood <- function(Assessment, figure = TRUE, save_figure = TRUE,
                               save_dir = getwd(), ...) {
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


#' Perform retrospective analysis
#'
#' Perform a retrospective analysis.
#'
#' @param Assessment An S4 object of Assessment output..
#' @param nyr The maximum number of years to remove for the retrospective analysis.
#' @param figure Indicates whether plots will be drawn.
#' @return A list with an array of model output and of model estimates from
#' the retrospective analysis.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @author Q. Huynh
#' @return Figures showing the time series of biomass and exploitation and parameter estimates
#' with successive number of years removed. Returns invisibly an array of model output and of model estimates.
#' @examples
#' \dontrun{
#' output <- DD_TMB(Data = DLMtool::Red_snapper)
#' get_retro <- retrospective(output, nyr = 5)
#' }
#' @export
retrospective <- function(Assessment, nyr = 5, figure = TRUE, save_figure = TRUE,
                          save_dir = getwd()) {
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


