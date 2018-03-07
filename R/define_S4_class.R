# Create "sdreport" class for output from TMB function sdreport()
setOldClass("sdreport")


#' Assessment Class
#'
#' An S4 class that contains objects from a function of class Assess.
#'
#' @name Assessment-class
#' @docType class
#'
#' @slot Model Name of the assessment model.
#' @slot MSY Maximum sustainable yield.
#' @slot UMSY Exploitation at maximum sustainable yield.
#' @slot FMSY Instantaneous fishing mortality rate at maximum sustainable yield.
#' @slot BMSY Biomass at maximum sustainable yield.
#' @slot B0 Biomass at virgin equilibrium.
#' @slot R0 Recruitment at virgin equilibrium.
#' @slot N0 Abundance at virgin equilibrium.
#' @slot SSB0 Spawning stock biomass at virgin equilibrium.
#' @slot h Steepness.
#' @slot U Time series of exploitation.
#' @slot U_UMSY Time series of relative exploitation.
#' @slot F Time series of fishing mortality.
#' @slot F_FMSY Time series of fishing mortality relative to MSY.
#' @slot B Time series of biomass.
#' @slot B_BMSY Time series of relative biomass relative to MSY.
#' @slot B_B0 Time series of depletion.
#' @slot SSB Time series of spawning stock biomass.
#' @slot SSB_SSBMSY Time series of spawning stock biomass relative to MSY.
#' @slot SSB_SSB0 Time series of spawning stock depletion.
#' @slot N Time series of population abundance.
#' @slot R Time series of recruitment.
#' @slot N_at_age Time series of numbers-at-age.
#' @slot C_at_age Time series of catch-at-age.
#' @slot selectivity Selectivity-at-age.
#' @slot Catch Predicted catch.
#' @slot Index Predicted Index.
#' @slot Random A vector of estimated random effects.
#' @slot Random_SE A vector of standard errors of the random effects.
#' @slot NLL Negative log-likelihood of the model (integrated across random effects).
#' @slot NLL_Catch Negative log-likelihood of the catch component.
#' @slot NLL_Index Negative log-likelihood of the index component.
#' @slot NLL_C_at_age Negative log-likelihood of the catch-at-age component.
#' @slot NLL_Random Marginal negative log-likelihood of the random effects.
#' @slot info A list containing the data and starting values of estimated parameters
#' for the assessment.
#' @slot obj A list with components returned from \code{\link[TMB]{MakeADFun}}.
#' @slot opt list with components from calling \code{\link[stats]{nlminb}} to \code{obj}.
#' @slot SD Parameter estimates and their standard errors, obtained from
#' \code{\link[TMB]{sdreport}}.
#' @slot TMB_report A list of model output reported from the TMB executable, i.e. \code{obj$report()}.
#' @slot dependencies A character string of data types used for the assessment.
#' @author Q. Huynh
#' @exportClass Assessment
setClass("Assessment", slots = c(Model = "character", MSY = "numeric",
                                 UMSY = "numeric", FMSY = "numeric", BMSY = "numeric",
                                 B0 = "numeric", R0 = "numeric", N0 = "numeric",
                                 SSB0 = "numeric", h = "numeric", U = "vector",
                                 U_UMSY = "numeric", F = "numeric", F_FMSY  = "numeric",
                                 B = "numeric", B_BMSY = "numeric", B_B0 = "numeric",
                                 SSB = "numeric", SSB_SSBMSY = "numeric", SSB_SSB0 = "numeric",
                                 N = "numeric", R = "numeric", N_at_age = "matrix",
                                 C_at_age = "matrix", selectivity = "numeric",
                                 Catch = "numeric", Index = "numeric",
                                 Random = "numeric", Random_SE = "numeric",
                                 NLL = "numeric", NLL_Catch = "numeric",
                                 NLL_Index = "numeric", NLL_C_at_age = "numeric",
                                 NLL_Random = "numeric",
                                 info = "list", obj = "list",
                                 opt = "list", SD = "sdreport", TMB_report = "list",
                                 dependencies = "character", Data = "Data"))


#' Summary of Assessment object
#'
#' @rdname summary-Assessment
#' @param object An object of class \linkS4class{Assessment}
#' @return A list of parameters
#' @examples
#' data(Red_snapper)
#' output <- DD_TMB(Red_snapper)
#' summary(output)
#' @exportMethod summary
setMethod("summary", signature(object = "Assessment"), function(object) {
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
#' data(Red_snapper)
#' output <- DD_TMB(Red_snapper)
#' plot(output, save_figure = FALSE)
#'
#' @exportMethod plot
setMethod("plot", signature(x = "Assessment"), function(x, save_figure = TRUE, save_dir = getwd()) {
  old.warning <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old.warning))

  f <- get(paste0("generate_plots_", x@Model))
  f(x, save_figure = save_figure, save_dir = save_dir)
})
