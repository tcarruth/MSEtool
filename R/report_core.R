

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
#' data(sim_snapper)
#' output <- DD_TMB(Data = sim_snapper)
#' pro <- profile_likelihood(output, UMSY = seq(0.05, 0.15, 0.01), MSY = seq(0.1, 3, 0.01))
#'
#' # Ensure your grid is of proper resolution.
#' # This grid is too coarse.
#' poor_profile <- profile_likelihood(output, UMSY = seq(0.05, 0.15, 0.1), MSY = seq(0.1, 3, 0.1),
#' save_figure = FALSE)
#' }
#' @export profile_likelihood
profile_likelihood <- function(Assessment, figure = TRUE, save_figure = TRUE,
                               save_dir = getwd(), ...) {
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
#' data(sim_snapper)
#' output <- DD_TMB(Data = sim_snapper)
#' get_retro <- retrospective(output, nyr = 5)
#' }
#' @export
retrospective <- function(Assessment, nyr = 5, figure = TRUE, save_figure = TRUE,
                          save_dir = getwd()) {
  f <- get(paste0('retrospective_', Assessment@Model))
  f(Assessment, nyr, figure = figure, save_figure = save_figure, save_dir = save_dir)
}


