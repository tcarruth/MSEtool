# Create "sdreport" class for output from TMB function sdreport()

#' @import DLMtool
#' @import graphics
#' @import stats
#' @import utils
#' @import methods
setOldClass("sdreport")
setOldClass("spictcls")
setClassUnion("sdreportAssess", members = c("character", "sdreport"))
setClassUnion("optAssess", members = c("list", "character"))
setClassUnion("listspict", members = c("list", "spictcls"))


#' Class-\code{Assessment}
#'
#' An S4 class that contains assessment output. Created from a function of class \code{Assess}.
#'
#' @name Assessment-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot Name Name of Data object.
#' @slot conv Logical. Whether the assessment model converged (defined by whether TMB returned a
#' positive-definite covariance matrix for the model).
#' @slot UMSY Estimate of exploitation at maximum sustainable yield.
#' @slot FMSY Estimate of instantaneous fishing mortality rate at maximum sustainable yield.
#' @slot MSY Estimate of maximum sustainable yield.
#' @slot BMSY Biomass at maximum sustainable yield.
#' @slot SSBMSY Spawning stock biomass at maximum sustainable yield.
#' @slot VBMSY Vulnerable biomass at maximum sustainable yield.
#' @slot B0 Biomass at unfished equilibrium.
#' @slot R0 Recruitment at unfished equilibrium.
#' @slot N0 Abundance at unfished equilibrium.
#' @slot SSB0 Spawning stock biomass at unfished equilibrium.
#' @slot VB0 Vulnerable biomass at unfished equilibrium.
#' @slot h Steepness.
#' @slot U Time series of exploitation.
#' @slot U_UMSY Time series of relative exploitation.
#' @slot FMort Time series of instantaneous fishing mortality.
#' @slot F_FMSY Time series of fishing mortality relative to MSY.
#' @slot B Time series of biomass.
#' @slot B_BMSY Time series of biomass relative to MSY.
#' @slot B_B0 Time series of depletion.
#' @slot SSB Time series of spawning stock biomass.
#' @slot SSB_SSBMSY Time series of spawning stock biomass relative to MSY.
#' @slot SSB_SSB0 Time series of spawning stock depletion.
#' @slot VB Time series of vulnerable biomass.
#' @slot VB_VBMSY Time series of vulnerable biomass relative to MSY.
#' @slot VB_VB0 Time series of vulnerable biomass depletion.
#' @slot R Time series of recruitment.
#' @slot N Time series of population abundance.
#' @slot N_at_age Time series of numbers-at-age matrix.
#' @slot Selectivity Selectivity-at-age matrix.
#' @slot Obs_Catch Observed catch.
#' @slot Obs_Index Observed index.
#' @slot Obs_C_at_age Observed catch-at-age matrix.
#' @slot Catch Predicted catch.
#' @slot Index Predicted index.
#' @slot C_at_age Predicted catch-at-age matrix.
#' @slot Dev A vector of estimated deviation parameters.
#' @slot Dev_type A description of the deviation parameters, e.g. "log recruitment deviations".
#' @slot NLL Negative log-likelihood. A vector for the total likelihood, integrated across random effects if applicable, components,
#' and penalty term (applied when \code{U > 0.975} in any year).
#' @slot SE_UMSY Standard error of UMSY estimate.
#' @slot SE_FMSY Standard error of FMSY estimate.
#' @slot SE_MSY Standard error of MSY estimate.
#' @slot SE_U_UMSY_final Standard error of U/UMSY in the terminal year.
#' @slot SE_F_FMSY_final Standard error of F/FMSY in the terminal year.
#' @slot SE_B_BMSY_final Standard error of B/BMSY in the terminal year.
#' @slot SE_B_B0_final Standard error of B/B0 in the terminal year.
#' @slot SE_SSB_SSBMSY_final Standard error of SSB/SSBMSY in the terminal year.
#' @slot SE_SSB_SSB0_final Standard error of SSB/SSB0 in the terminal year.
#' @slot SE_VB_VBMSY_final Standard error of VB/VBMSY in the terminal year.
#' @slot SE_VB_VB0_final Standard error of VB/VB0 in the terminal year.
#' @slot SE_Dev A vector of standard errors of the deviation parameters.
#' @slot info A list containing the data and starting values of estimated parameters
#' for the assessment.
#' @slot obj A list with components returned from \code{\link[TMB]{MakeADFun}}.
#' @slot opt A list with components from calling \code{\link[stats]{nlminb}} to \code{obj}.
#' @slot SD A list (class sdreport) with parameter estimates and their standard errors, obtained from
#' \code{\link[TMB]{sdreport}}.
#' @slot TMB_report A list of model output reported from the TMB executable, i.e. \code{obj$report()}, and derived quantities (e.g. MSY).
#' @slot dependencies A character string of data types required for the assessment.
#' @examples
#' \donttest{
#' output <- DD_TMB(Data = DLMtool::Red_snapper)
#' class(output)
#' }
#' @seealso \link{plot.Assessment} \link{summary.Assessment} \link{retrospective} \link{profile_likelihood} \link{make_MP}
#' @author Q. Huynh
#' @export Assessment
#' @exportClass Assessment
Assessment <- setClass("Assessment",
                       slots = c(Model = "character", Name = "character", conv = "logical", UMSY = "numeric", FMSY = "numeric",
                                 MSY = "numeric", BMSY = "numeric", SSBMSY = "numeric", VBMSY = "numeric",
                                 B0 = "numeric", R0 = "numeric", N0 = "numeric", SSB0 = "numeric", VB0 = "numeric",
                                 h = "numeric", U = "numeric", U_UMSY = "numeric", FMort = "numeric", F_FMSY  = "numeric",
                                 B = "numeric", B_BMSY = "numeric", B_B0 = "numeric",
                                 SSB = "numeric", SSB_SSBMSY = "numeric", SSB_SSB0 = "numeric", VB = "numeric",
                                 VB_VBMSY = "numeric", VB_VB0 = "numeric",
                                 R = "numeric", N = "numeric", N_at_age = "matrix",
                                 Selectivity = "matrix", Obs_Catch = "numeric", Obs_Index = "numeric",
                                 Obs_C_at_age = "matrix", Catch = "numeric", Index = "numeric",
                                 C_at_age = "matrix", Dev = "numeric", Dev_type = "character",
                                 NLL = "numeric", SE_UMSY = "numeric", SE_FMSY = "numeric", SE_MSY = "numeric",
                                 SE_U_UMSY_final = "numeric", SE_F_FMSY_final = "numeric",
                                 SE_B_BMSY_final = "numeric", SE_B_B0_final = "numeric",
                                 SE_SSB_SSBMSY_final = "numeric", SE_SSB_SSB0_final = "numeric",
                                 SE_VB_VBMSY_final = "numeric", SE_VB_VB0_final = "numeric",
                                 SE_Dev = "numeric", info = "listspict", obj = "list", opt = "optAssess", SD = "sdreportAssess",
                                 TMB_report = "list", dependencies = "character"))

#' @name summary.Assessment
#' @title Summary of Assessment object
#' @aliases summary,Assessment-method
#' @description Returns a summary of parameter estimates and output from an \linkS4class{Assessment} object.
#' @param object An object of class \linkS4class{Assessment}
#' @return A list of parameters.
#' @examples
#' output <- DD_TMB(Data = DLMtool::Simulation_1)
#' summary(output)
#' @exportMethod summary
setMethod("summary", signature(object = "Assessment"), function(object) {
  if(is.character(object@opt) || is.character(object@SD)) warning("Did model converge?")
  f <- get(paste0("summary_", object@Model))
  f(object)
})

#' @name plot.Assessment
#' @aliases plot,Assessment,missing-method plot,Assessment,retro
#' @title Plot Assessment object
#' @description Produces HTML file (via markdown) figures of parameter estimates and output from an \linkS4class{Assessment} object.
#'
#' @param x An object of class \linkS4class{Assessment}.
#' @param y An object of class \linkS4class{retro}.
#' @param filename Character string for the name of the markdown and HTML files.
#' @param dir The directory in which the markdown and HTML files will be saved.
#' @param ret_yr If greater than zero, then a retrospective analysis will be performed and results will be reported. The integer here corresponds
#' to the number of peels (the maximum number of terminal years for which the data are removed).
#' @param open_file Logical, whether the HTML document is opened after it is rendered.
#' @param quiet Logical, whether to silence the markdown rendering function.
#' @param ... Other arguments to pass to \link[rmarkdown]{render}.
#' @return Returns invisibly the output from \link[rmarkdown]{render}.
#' @examples
#' \donttest{
#' output <- DD_TMB(Data = Simulation_1)
#' plot(output)
#' }
#' @seealso \link{retrospective}
#' @exportMethod
setMethod("plot", signature(x = "Assessment", y = "missing"),
          function(x, filename = paste0("report_", x@Model), dir = tempdir(), ret_yr = 0L, open_file = TRUE, quiet = TRUE, ...) {
            # Generating retrospective
            if(is.numeric(ret_yr) && ret_yr > 0) {
              ret_yr <- as.integer(ret_yr)
              message("Running retrospective...")
              ret <- retrospective(x, nyr = ret_yr, figure = FALSE)
            } else ret <- NULL
            report(x, ret, filename = filename, dir = dir, open_file = open_file, quiet = quiet, ...)
          })

#' @name retrospective
#' @title Retrospective analysis of assessment models
#'
#' @description Perform a retrospective analysis, successive removals of most recent years of data to evaluate resulting
#' parameter estimates.
#'
#' @param x An S4 object of class \linkS4class{Assessment}.
#' @param nyr The maximum number of years to remove for the retrospective analysis.
#' @param figure Indicates whether plots will be drawn.
#' @return A list with an array of model output and of model estimates from
#' the retrospective analysis.
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
setGeneric("retrospective", function(x, ...) standardGeneric("retrospective"))

#' @rdname retrospective
#' @aliases retrospective,Assessment-method
#' @exportMethod
setMethod("retrospective", signature(x = "Assessment"),
          function(x, nyr = 5, figure = TRUE) {
            if(figure) {
              old.warning <- options()$warn
              options(warn = -1)
              on.exit(options(warn = old.warning))

              old_par <- par(no.readonly = TRUE)
              on.exit(par(old_par), add = TRUE)
            }

            f <- get(paste0('retrospective_', x@Model))
            f(x, nyr, figure)
          })


#' Class-\code{retro}
#'
#' An S4 class that contains output from \link{retrospective}.
#'
#' @name retro-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot Name Name of Data object.
#' @slot TS_var Character vector of time series variables, e.g. recruitment, biomass, from the assessment.
#' @slot TS An array of time series assessment output of dimension, indexed by: peel (the number of terminal years removed from the base assessment),
#' years, and variables (corresponding to `TS_var`).
#' @slot Est_var Character vector of estimated paramters, e.g. R0, steeppness, in the assessment.
#' @slot Est An array for estimated parameters of dimension, indexed by: peel, variables (corresponding to `Est_var`), and
#' value (length 2 for estimate and standard error).
#' @seealso \link{plot.retro} \link{summary.retro} \link{plot.Assessment}
#' @author Q. Huynh
#' @export retro
#' @exportClass retro
retro <- setClass("retro", slots = c(Model = "character", Name = "character", TS_var = "character",
                                     TS = "array", Est_var = "character", Est = "array"))

#' @rdname plot.Assessment
#' @aliases plot,Assessment,retro-method
#' @exportMethod
setMethod("plot", signature(x = "Assessment", y = "retro"),
          function(x, y, filename = paste0("report_", x@Model), dir = tempdir(), open_file = TRUE, quiet = TRUE, ...) {
            report(x, y, filename = filename, dir = dir, open_file = open_file, quiet = quiet, ...)
          })

#' @rdname retrospective
#' @aliases plot.retro plot,retro,missing-method
#' @exportMethod
setMethod("plot", signature(x = "retro", y = "missing"),
          function(x, color = NULL) {
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))

            retro <- x
            if(is.null(color)) color <- rich.colors(dim(retro@TS)[1])
            Year_matrix <- matrix(as.numeric(dimnames(retro@TS)$Year), ncol = length(color), nrow = dim(retro@TS)[2], byrow = FALSE)
            xlim <- range(as.numeric(dimnames(retro@TS)$Year))
            nyr_label <- dimnames(retro@TS)$Peel

            for(i in 1:length(retro@TS_var)) {
              matrix_to_plot <- t(retro@TS[, , i])
              ylim <- c(0, 1.1 * max(matrix_to_plot, na.rm = TRUE))
              ylab <- attr(retro, "TS_lab")[i]

              plot(NULL, NULL, xlim = xlim, ylim = ylim, xlab = "Year", ylab = ylab)
              abline(h = 0, col = "grey")
              if(grepl("MSY", as.character(ylab))) abline(h = 1, lty = 3)

              matlines(Year_matrix, matrix_to_plot, col = color, lty = 1)

              legend("topleft", legend = nyr_label, lwd = 1, col = color, bty = "n", title = "Years removed:")
            }
            invisible()
          })


setMethod("show", signature(object = "retro"), function(object) print(calculate_Mohn_rho(object@TS, ts_lab = attr(object, "TS_lab"))))

#' @rdname retrospective
#' @aliases summary.retro summary,retro-method
#' @exportMethod
setMethod("summary", signature(object = "retro"), function(object) calculate_Mohn_rho(object@TS, ts_lab = attr(object, "TS_lab")))


if(getRversion() >= "2.15.1") {
  # Define global variables for Assessment objects
  utils::globalVariables(slotNames("Assessment"))

  utils::globalVariables("plot.dir")

  # For Awatea2OM - Quang assumes these variables are loaded in from .rda files
  utils::globalVariables(c("Bmcmc", "currentMCMC", "currentRes"))
}


