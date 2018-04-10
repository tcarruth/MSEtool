# ===================================================================
# === Harvest control rules =========================================
# ===================================================================

#' A harvest control rule to fish at maximum sustainable yield.
#'
#' A simple control rule that specifies the total allowable catch (TAC) to be the
#' product of current vulnerable biomass (VB_final) and FMSY (or UMSY).
#'
#' @param Assessment An object of class \linkS4class{Assessment} with estimates of
#' FMSY or UMSY and vulnerable biomass in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param ... Miscellaneous arguments (not currently used).
#' @return An object of class \linkS4class{Rec} with the TAC recommendation.
#' @author Q. Huynh
#' @references
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#' @seealso \link{make_MP}
#' @examples
#' # create an MP to run in closed-loop MSE
#' DD_MSY <- make_MP(DD_TMB, HCR_MSY)
#' class(DD_MSY)
#' \dontrun{
#' myOM <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "DD_MSY"))
#' }
#' @export
HCR_MSY <- function(Assessment, reps = 1, ...) {
  TAC <- calculate_TAC(Assessment, reps)
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR_MSY) <- "HCR"



#' Ramped harvest control rules (40-10 and 60-20)
#'
#' A output control rule with a ramp that reduces the TAC recommendation
#' (FMSY * Vulnerable biomass) when spawning biomass depletion (SSB/SSB0) is less than
#' the target reference point (TRP). The TAC reduction is linearly reduced to 0
#' when depletion is less than the limit reference point (LRP). For example,
#' TRP is 0.4 and LRP is 0.1 in the 40-10 control rule.
#'
#' @param Assessment An object of class \linkS4class{Assessment} with estimates of
#' FMSY or UMSY, vulnerable biomass, and spawning biomass depletion in terminal year.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param ... Miscellaneous arguments (not currently used).
#' @note Function will use slot B_B0 of the Asssessment object (e.g. for surplus
#' production models \link{SP} and \link{SP_SS}) if SSB_SSB0 is unavailable.
#' @return An object of class \linkS4class{Rec} with the TAC recommendation.
#' @author Q. Huynh
#' @references
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#' @seealso \link{HCR60_20} \link{HCRlin} \link{make_MP}
#' @examples
#' # 40-10 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.1, 0.4), xlab = "Estimated SSB/SSB0",
#' ylab = "TAC adjustment factor", main = "40-10 harvest control rule",
#' type = "l", col = "blue")
#' abline(v = c(0.1, 0.4), col = "red", lty = 2)
#'
#' # create an MP to run in closed-loop MSE
#' DD_40_10 <- make_MP(DD_TMB, HCR40_10)
#'
#' \dontrun{
#' myOM <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "DD_40_10"))
#' }
#'
#' # 60-20 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.2, 0.6), xlab = "Estimated SSB/SSB0",
#' ylab = "TAC adjustment", main = "60-20 harvest control rule", type = 'l',
#' col = 'blue')
#' abline(v = c(0.2, 0.6), col = "red", lty = 2)
#'
#' # create an MP to run in closed-loop MSE
#' DD_60_20 <- make_MP(DD_TMB, HCR60_20)
#'
#' @describeIn HCR40_10 Common U.S. west coast control rule.
#' @export
HCR40_10 <- function(Assessment, reps = 1, ...) {
  TAC <- calculate_TAC(Assessment, reps)
  if(!all(is.na(TAC))) {
    if(length(Assessment@SSB_SSB0) > 0) {
      depletion <- Assessment@SSB_SSB0[length(Assessment@SSB_SSB0)]
    } else
      if(length(Assessment@B_B0) > 0) {
        depletion <- Assessment@B_B0[length(Assessment@B_B0)]
      } else depletion <- NA
    ramp <- HCRlin(depletion, 0.1, 0.4)
    TAC <- ramp * TAC
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR40_10) <- "HCR"



#' @describeIn HCR40_10 More conservative than 40-10
#' @export
HCR60_20 <- function(Assessment, reps = 1, ...) {
  TAC <- calculate_TAC(Assessment, reps)
  if(!all(is.na(TAC))) {
    if(length(Assessment@SSB_SSB0) > 0) {
      depletion <- Assessment@SSB_SSB0[length(Assessment@SSB_SSB0)]
    } else
      if(length(Assessment@B_B0) > 0) {
        depletion <- Assessment@B_B0[length(Assessment@B_B0)]
      } else depletion <- NA
    ramp <- HCRlin(depletion, 0.2, 0.6)
    TAC <- ramp * TAC
  }
  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR60_20) <- "HCR"



#' Generic linear harvest control rule based on biomass
#'
#' A general function that adjusts the TAC by a linear ramp based on estimated biomass.
#'
#' @param Brel Improper fraction: An estimate of biomass (either absolute
#' or relative, e.g. B/BMSY or B/B0).
#' @param LRP Improper fraction: the Limit Reference Point, the biomass
#' below which the adjustment is zero (no fishing). Same units as \code{Brel}.
#' @param TRP Improper fraction: the Target Reference Point, the biomass
#' above which the adjustment is 1 (no adjustment). Same units as \code{Brel}.
#' @return a TAC or TAE adjustement factor.
#' @author T. Carruthers
#' @references
#' Deroba, J.J. and Bence, J.R. 2008. A review of harvest policies: Understanding relative
#' performance of control rules. Fisheries Research 94:210-223.
#'
#' Restrepo, V.R. and Power, J.E. 1999. Precautionary control rules in US fisheries
#' management: specification and performance. ICES Journal of Marine Science 56:846-852.
#' @export HCRlin
#' @examples
#' #40-10 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.1, 0.4), xlab = "Estimated B/B0", ylab = "TAC adjustement factor",
#' main = "A 40-10 harvest control rule", type = 'l', col = 'blue')
#' abline(v = c(0.1,0.4), col = 'red', lty = 2)
HCRlin <- function(Brel, LRP, TRP){
  adj <- rep(1, length(Brel))
  adj[Brel <= LRP] <- 0
  cond <- Brel>LRP & Brel<TRP
  adj[cond] <- (Brel[cond] - LRP)/(TRP - LRP)
  adj
}


#' A Harvest Control Rule using B/BMSY and F/FMSY to adjust TAC or TAE.
#'
#' @param Brel improper fraction: an estimate of Biomass relative to BMSY
#' @param Frel improper fraction: an estimate of Fishing mortality rate relative to FMSY
#' @param Bpow non-negative real number: controls the shape of the biomass adjustment, when zero there is no adjustment
#' @param Bgrad non-negative real number: controls the gradient of the biomass adjustment
#' @param Fpow non-negative real number: controls the adjustment speed relative to F/FMSY. When set to 1, next recommendation is FMSY. When less than 1 next recommendation is between current F and FMSY.
#' @param Fgrad improper fraction: target Fishing rate relative to FMSY
#' @return a TAC or TAE adjustement factor.
#' @author T. Carruthers
#' @references Made up for this package
#' @export
#' @examples
#' res <- 100
#' Frel <- seq(1/2, 2, length.out = res)
#' Brel <- seq(0.05, 2, length.out=res)
#' adj <- array(HCR_FB(Brel[rep(1:res, res)], Frel[rep(1:res, each = res)],
#'                     Bpow = 2, Bgrad = 1, Fpow = 1, Fgrad = 0.75), c(res, res))
#' contour(Brel, Frel, adj, nlevels = 20, xlab = "B/BMSY", ylab = "F/FMSY",
#'         main = "FBsurface TAC adjustment factor")
#' abline(h = 1, col = 'red', lty = 2)
#' abline(v = 1, col = 'red', lty = 2)
#' legend('topright', c("Bpow = 2", "Bgrad = 1", "Fpow = 1", "Fgrad = 0.75"), text.col = 'blue')
HCR_FB<-function(Brel,Frel,Bpow=2,Bgrad=1,Fpow=1,Fgrad=1){
  Fresp <- exp(log(1/Frel)*Fpow)
  Bresp <- exp(powdif(Brel-1,Bpow,Bgrad))
  Fgrad*Bresp*Fresp
}


#' Power difference function.
#'
#' @param x Real number: the absolute difference between two numbers
#' @param z Real number: the exponent the difference will be raised to
#' @param g Real number: the gradient in the exponential difference
#' @return a positive real number
#' @export
#' @examples
#' powdif(-2,3,1)
powdif<-function(x,z,g){
  x2<-(g*abs(x))^z
  x2[x<0]<-(-x2[x<0])
  x2
}

#' Calculate TAC from Assessment object
#'
#' A function to calculate the total allowable catch TAC as the product of
#' either UMSY or FMSY and the vulnerable biomass in terminal year.
#'
#' @param Assessment An Assessment object with estimates of UMSY or FMSY and
#' terminal year vulnerable biomass.
#' @param reps The number of stochastic draws of UMSY or FMSY (Set equal to one
#' for deterministic calculations).
#' @return A vector of length \code{reps} of stochastic samples of TAC recommendation.
#' @seealso \link{HCR_MSY} \link{HCR40_10} \link{HCR60_20}
#' @export
calculate_TAC <- function(Assessment, reps) {
  has_UMSY <- has_FMSY <- has_VB <- FALSE
  if(length(Assessment@UMSY) > 0) has_UMSY <- TRUE
  if(length(Assessment@FMSY) > 0) has_FMSY <- TRUE
  if(length(Assessment@VB) > 0) has_VB <- TRUE

  if(has_VB && (has_UMSY || has_FMSY)) {
    VB_current <- Assessment@VB[length(Assessment@VB)]
    if(has_UMSY) {
      if(length(Assessment@SE_UMSY) > 0) {
        SE_UMSY <- Assessment@SE_UMSY
      } else SE_UMSY <- 1e-8
      UMSY_vector <- trlnorm(reps, Assessment@UMSY, SE_UMSY)
      TAC <- UMSY_vector * VB_current
    }
    if(has_FMSY) {
      if(length(Assessment@SE_FMSY) > 0) {
        SE_FMSY <- Assessment@SE_FMSY
      } else SE_FMSY <- 1e-8
      FMSY_vector <- trlnorm(reps, Assessment@FMSY, SE_FMSY)
      TAC <- FMSY_vector * VB_current
    }
  } else {
    TAC <- rep(NA, reps) # Missing estimates for HCR.
  }
  return(TAC)
}
