
#' Generate assessment report
#'
#' Returns a list of assessment output and, optionally, a multitude of plots
#' (e.g., residuals, predicted time series).
#'
#' @param Assessment An S4 object of assessment output generated from an MP.
#' @param figure Indicates whether plots will be drawn.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @return A list of data and parameter inputs, predicted time series and derived quantities,
#' and model estimates and standard deviations.
#' Optionally, corresponding plots will also be drawn.
#' @export generate_report
#' @details List dataframes and figures here.
#' @author Q. Huynh
#' @examples
#' dir <- tempdir()
#' setwd(dir)
#'
#' output <- DD_TMB(1, Red_snapper, Assessment = TRUE)
#' RS_report <- generate_report(output, figure = FALSE, save_dir = dir)
generate_report <- function(Assessment, figure = TRUE, save_figure = TRUE,
                            save_dir = getwd()) {
  old.stringsAsFactors <- options()$stringsAsFactors
  on.exit(options(stringsAsFactors = old.stringsAsFactors))
  options(stringsAsFactors = FALSE)

  f <- get(paste0('generate_report_', Assessment@MP))
  f(Assessment, figure = figure, save_figure = save_figure, save_dir = save_dir)
}

#' Perform profile likelihood
#'
#' Performs a profile likelihood.
#'
#' @param Assessment An S4 object of assessment output generated from an MP.
#' @param figure Indicates whether a figure will be plotted.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @param ... A sequence of values of the parameter(s) for the profile.
#' See details for name of arguments to be passed on.
#' @details For \code{DD_TMB}, sequence of values for \code{UMSY} and \code{MSY} are needed.
#' @author Q. Huynh
#' @return A data frame of negative log-likelihood values from profile and, optionally,
#' a figure.
#' @export profile_likelihood
profile_likelihood <- function(Assessment, figure = TRUE, save_figure = TRUE,
                               save_dir = getwd(), ...) {
  f <- get(paste0('profile_likelihood_', Assessment@MP))
  f(Assessment, figure = figure, save_figure = save_figure, save_dir = save_dir, ...)
}


#' Perform retrospective analysis
#'
#' Perform a retrospective analysis.
#'
#' @param Assessment An S4 object of assessment output generated from an MP.
#' @param nyr The maximum number of years to remove for the retrospective analysis.
#' @param figure Indicates whether plots will be drawn.
#' @return A list with an array of model output and of model estimates from
#' the retrospective analysis.
#' @param save_figure Indicates whether figures will be saved to directory.
#' @param save_dir The directory to which figures will be saved. By default: \code{getwd()}
#' @export retrospective
#' @author Q. Huynh
#' @examples
#' output <- DD_TMB(1, Red_snapper, Assessment = TRUE)
#' RS_retro <- retrospective(output, nyr = 5)
retrospective <- function(Assessment, nyr = 5, figure = TRUE, save_figure = TRUE,
                          save_dir = getwd()) {
  f <- get(paste0('retrospective_', Assessment@MP))
  f(Assessment, nyr, figure = figure, save_figure = save_figure, save_dir = save_dir)
}


