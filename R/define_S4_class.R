# Create "sdreport" class for output from TMB function sdreport()

#' @import DLMtool
#' @import graphics
#' @import stats
#' @import utils
#' @import methods
setOldClass("sdreport")
setClassUnion("sdreportAssess", members = c("character", "sdreport"))
setClassUnion("optAssess", members = c("list", "character"))



#' Class-\code{Assessment}
#'
#' An S4 class that contains assessment output. Created from a function of class \code{Assess}.
#'
#' @name Assessment-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot UMSY Estimate of exploitation at maximum sustainable yield.
#' @slot FMSY Estimate of instantaneous fishing mortality rate at maximum sustainable yield.
#' @slot MSY Estimate of maximum sustainable yield.
#' @slot BMSY Biomass at maximum sustainable yield.
#' @slot SSBMSY Spawning stock biomass at maximum sustainable yield.
#' @slot VBMSY Vulnerable biomass at maximum sustainable yield.
#' @slot B0 Biomass at virgin equilibrium.
#' @slot R0 Recruitment at virgin equilibrium.
#' @slot N0 Abundance at virgin equilibrium.
#' @slot SSB0 Spawning stock biomass at virgin equilibrium.
#' @slot VB0 Vulnerable biomass at virgin equilibrium.
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
#' @slot NLL Negative log-likelihood (total [integrated across random effects] and components).
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
#' @slot opt list with components from calling \code{\link[stats]{nlminb}} to \code{obj}.
#' @slot SD Parameter estimates and their standard errors, obtained from
#' \code{\link[TMB]{sdreport}}.
#' @slot TMB_report A list of model output reported from the TMB executable, i.e. \code{obj$report()}.
#' @slot dependencies A character string of data types used for the assessment.
#' @slot Data An object of class \linkS4class{Data} that was used to perform the assessment.
#' @examples
#' \dontrun{
#' data(sim_snapper)
#' output <- DD_TMB(1, sim_snapper)
#' class(output)
#' str(output)
#' }
#' @author Q. Huynh
#' @export
#' @exportClass Assessment
Assessment <- setClass("Assessment",
 slots = c(Model = "character", UMSY = "numeric", FMSY = "numeric",
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
 SE_Dev = "numeric", info = "list", obj = "list", opt = "optAssess", SD = "sdreportAssess",
 TMB_report = "list", dependencies = "character", Data = "Data"))


#' Summary of Assessment object
#'
#' @rdname summary-Assessment
#' @param object An object of class \linkS4class{Assessment}
#' @return A list of parameters
#' @examples
#' \dontrun{
#' output <- DD_TMB(Data = DLMtool::Simulation_1)
#' summary(output)
#' }
#' @exportMethod summary
setMethod("summary", signature(object = "Assessment"), function(object) {
  if(is.character(object@opt) || is.character(object@SD)) warning("Did model converge? Check slots obj, opt, and SD.")
  f <- get(paste0("summary_", object@Model))
  f(object)
})

#' Plot Assessment object
#'
#' @rdname plot-Assessment
#' @param x An object of class \linkS4class{Assessment}
#' @param save_figure Indicates whether figures will be saved to disk. A corresponding
#' html report will be produced.
#' @param save_dir The directory (by default, the current working directory) in which a
#' sub-directory will be created to save figures.
#' @examples
#' \dontrun{
#' output <- DD_TMB(Data = Simulation_1)
#' plot(output, save_figure = FALSE)
#' }
#' @exportMethod plot
setMethod("plot", signature(x = "Assessment"), function(x, save_figure = TRUE, save_dir = getwd()) {
  if(is.character(x@opt) || is.character(x@SD)) warning("Did model converge? Check slots obj, opt, and SD.")
  old.warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old.warning))
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  f <- get(paste0("generate_plots_", x@Model))
  f(x, save_figure = save_figure, save_dir = save_dir)
})


if(getRversion() >= "2.15.1") {
  # Define global variables for Assessment objects
  utils::globalVariables(slotNames("Assessment"))

  utils::globalVariables("plot.dir")

  # For Awatea2OM - Quang assumes these variables are loaded in from .rda files
  utils::globalVariables(c("Bmcmc", "currentMCMC", "currentRes"))
}


