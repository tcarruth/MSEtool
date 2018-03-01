#' Assessment Class
#'
#'
#' An S4 class that contains objects from an assessment called in an MP.
#'
#' @name Assessment-class
#' @docType class
#'
#' @slot MP Name of the MP.
#' @slot Rec An object of class \code{\link[DLMtool]\linkS4class{Rec}} which contains
#' the management recommendation from the MP.
#' @slot info A data frame containing the data and starting value of estimated parameters
#' (excluding nuisance parameters with analytic solutions) for the assessment.
#' @slot obj A list with components returned from \code{\link[TMB]{MakeADFun}}.
#' @slot opt list with components from calling \code{\link[stats]{nlminb}} to \code{obj}.
#' @slot SD Parameter estimates and their standard errors, obtained from
#' \code{\link[TMB]{sdreport}}.
#' @slot report A list of model output reported from TMB executable, i.e. \code{obj$report()}.
#' @slot Data An object of class \code{\link[DLMtool]\linkS4class{Data}} for
#' implementing the MP.
#' @slot dependencies A character string of data types used for the MP.
setClass("Assessment", representation(MP = "character", Rec = "Rec", info = "list", obj = "list",
                                      opt = "list", SD = "list", report = "list", Data = "Data",
                                      dependencies = "character"))
