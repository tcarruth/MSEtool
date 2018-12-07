# MOM object

# ---- MOM Class ----
#' Class \code{'MOM'}
#'
#' An object containing all the parameters needed to control a multi-stock, multi-fleet MSE which can
#' be build from component Stock, Fleet, Obs, and Imp objects.
#'
#' Almost all of these inputs are a vector of length 2 which describes the upper and lower
#' bounds of a uniform distribution from which to sample the parameter.
#'
#'
#' @name MOM-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new('MOM', Stock_list, Fleet_list, Obs_list, Imp_list)}.

#' @slot Name Name of the operating model
#' @slot Agency Name of the agency responsible for the management of the fishery. Character string
#' @slot Region Name of the general geographic region of the fishery. Character string
#' @slot Sponsor Name of the organization who sponsored the OM. Character string
#' @slot Latitude Latitude (decimal degrees). Negative values represent the South of the Equator. Numeric. Single value
#' @slot Longitude Longitude (decimal degrees). Negative values represent the West of the Prime Meridian. Numeric. Single value

#' @slot nsim The number of simulations
#' @slot proyears The number of projected years
#' @slot interval The assessment interval - how often would you like to update the management system?
#' @slot pstar The percentile of the sample of the management recommendation for each method
#' @slot maxF Maximum instantaneous fishing mortality rate that may be simulated for any given age class
#' @slot reps Number of samples of the management recommendation for each method. Note that when this is set to 1, the mean value of
#' the data inputs is used.
#' @slot Stocks_cpars A hierarcical list nstock long of custom parameters. Time series are a matrix nsim rows by nyears columns. Single parameters are a vector nsim long
#' @slot Fleets_cpars A hierarcical list nfleet long of custom parameters. Time series are a matrix nsim rows by nyears columns. Single parameters are a vector nsim long
#' @slot seed A random seed to ensure users can reproduce results exactly
#' @slot Source A reference to a website or article from which parameters were taken to define the operating model
#' @slot Stocks List of stock objects
#' @slot Fleets List of Fleet objects
#' @slot Obs    Hierarchical List of Observation model objects Level 1 is stock, level 2 is fleet
#' @slot Imps   Hierarchical List of Implementation model objects Level 1 is stock, level 2 is fleet
#' @slot CatchFrac A list nstock long, of vectors representing the fraction of current catches of the various fleets to each stock (each vector is nfleet long and sums to 1 for each stock)
#' @slot Allocation  A list nstock long, of vector allocations by fleet (default is NULL and allocation is according to last historical year). Note this
#' over-ridden if an MP of class 'MP_F" is supplied that is a multi-fleet MP.
#' @author T. Carruthers and A. Hordyk
#' @export
#' @keywords classes
#'
setClass("MOM", representation(Name = "character", Agency="character",
                              Region="character", Sponsor="character",
                              Latitude="numeric", Longitude="numeric",
                              nsim="numeric", proyears="numeric",
                              interval='numeric', pstar='numeric', maxF='numeric', reps='numeric',
                              cpars="list", seed="numeric", Source="character",Stocks='list',Fleets='list',Obs='list',Imps='list',
                              CatchFrac='list',Allocation='list'))


# initialize MOM
setMethod("initialize", "MOM", function(.Object, Stocks=NULL, Fleets=NULL, Obs=NULL, Imps=NULL, CatchFrac=NULL, Stocks_cpars=NULL,Fleets_cpars=NULL,
                                       interval=4, pstar=0.5, maxF=0.8, reps=1, nsim=48, proyears=50, Source=NULL, Allocation=NULL) {

  if (is.null(Stocks)|is.null(Fleets)|is.null(Obs)|is.null(Imps)|is.null(CatchFrac)) {
    message("A specified list of objects is required for each of the following arguments: Stocks, Fleets, Obs, Imps, CatchFrac. Returning a blank MOM object")
    .Object@seed <- 1
    return(.Object)
  }

  # helper function to expand Fleets by stock

  for(i in 1:length(Stocks))if (class(Stocks[[i]]) != "Stock")  stop(paste("Could not build operating model:", deparse(substitute(Stocks[[i]])), "not of class Stock"))
  for(i in 1:length(Fleets))for(j in 1:length(Fleets[[i]]))if (class(Fleets[[i]][[j]]) != "Fleet")  stop(paste("Could not build operating model:", deparse(substitute(Fleets[[i]][[j]])), "not of class Fleet"))
  for(i in 1:length(Obs))for(j in 1:length(Obs[[i]]))  if (class(Obs[[i]][[j]]) != "Obs") stop(paste("Could not build operating model:", deparse(substitute(Obs[[i]][[j]])), "not of class Obs"))
  for(i in 1:length(Imps))for(j in 1:length(Imps[[i]])) if (class(Imps[[i]][[j]]) != "Imp") stop(paste("Could not build operating model:", deparse(substitute(Imp)), "not of class Imp"))
  .Object@Name <- paste(c("MultiOM:",SIL(Stocks,"Name"),SIL(Fleets,"Name")),collapse=" - ")

  if( length(unique(SIL(Fleets,"nyears")))>1)stop("Fleet objects must have the same historical duration (nyears)")

   # Now copy the values for stock, fleet and observation slots to same
  .Object@Stocks<-Stocks
  .Object@Fleets<-Fleets
  .Object@cpars<-cpars
  .Object@Obs<-Obs
  .Object@Imps<-Imps

  if(is.null(Source)){
    .Object@Source <- Stocks[[1]]@Source
  }else{
    .Object@Source <- Source
  }

  .Object@nsim <- nsim
  .Object@proyears <- proyears
  .Object@interval <- interval
  .Object@pstar <- pstar
  .Object@maxF <- maxF
  .Object@reps <- reps
  .Object@Allocation <- new('list')
  .Object@CatchFrac=CatchFrac
  .Object@seed=1
  .Object

})
