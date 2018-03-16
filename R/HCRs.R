# ===================================================================
# === Harvest control rules =========================================
# ===================================================================

#' A harvest control rule to fish at maximum sustainable yield.
#'
#' A simple control rule that specifies the total allowable catch (TAC) to be the
#' product of current biomass and FMSY (or UMSY).
#'
#' @param Assessment An object of class \linkS4class{Assessment}.
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
  has_UMSY <- has_FMSY <- has_B <- FALSE
  if(!is.na(Assessment@UMSY[1])) has_UMSY <- TRUE
  if(!is.na(Assessment@FMSY[1])) has_FMSY <- TRUE
  if(length(Assessment@B) > 0) has_B <- TRUE

  if(has_B && (has_UMSY || has_FMSY)) {
    B_current <- Assessment@B[length(Assessment@B)]
    if(has_UMSY) {
      SE_UMSY <- ifelse(is.na(Assessment@SE_UMSY[1]), 1e-8, Assessment@SE_UMSY)
      UMSY_vector <- trlnorm(reps, Assessment@UMSY, SE_UMSY)
      TAC <- UMSY_vector * B_current
    }
    if(has_FMSY) {
      SE_FMSY <- ifelse(is.na(Assessment@SE_FMSY[1]), 1e-8, Assessment@SE_FMSY)
      FMSY_vector <- trlnorm(reps, Assessment@FMSY, SE_FMSY)
      TAC <- FMSY_vector * B_current
    }
  } else {
    TAC <- rep(NA, reps) # Missing estimates for HCR. Previous TAC recommendation will be used.
  }

  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR_MSY) <- "HCR"



#' A 40-10 harvest control rule.
#'
#' A output control rule with a ramp that reduces the
#' TAC recommendation (FMSY * Biomass) when biomass depletion (B/B0) is less than 40%.
#' The TAC reduction is linearly reduced to 0 when depletion is less than 10%.
#'
#' @param Assessment An object of class \linkS4class{Assessment}.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param ... Miscellaneous arguments (not currently used).
#' @return An object of class \linkS4class{Rec} with the TAC recommendation.
#' @author Q. Huynh
#' @references
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#' @seealso \link{HCR60_20} \link{HCRlin} \link{make_MP}
#' @examples
#' # 40-10 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.1, 0.4), xlab = "Estimated B/B0", ylab = "TAC adjustement factor",
#' main = "A 40-10 harvest control rule", type = "l", col = "blue")
#' abline(v = c(0.1,0.4), col = "red", lty = 2)
#'
#' # create an MP to run in closed-loop MSE
#' DD_40_10 <- make_MP(DD_TMB, HCR40_10)
#' \dontrun{
#' myOM <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "DD_40_10"))
#' }
#' @export
HCR40_10 <- function(x = 1, Assessment, reps = 1, ...) {
  has_UMSY <- has_FMSY <- has_B <- has_B_B0 <- FALSE
  if(!is.na(Assessment@UMSY[1])) has_UMSY <- TRUE
  if(!is.na(Assessment@FMSY[1])) has_FMSY <- TRUE
  if(length(Assessment@B) > 0) has_B <- TRUE
  if(length(Assessment@B_B0) > 0) has_B_B0 <- TRUE

  if((has_UMSY || has_FMSY) && has_B && has_B_B0) {
    B_B0 <- Assessment@B_B0[length(Assessment@B_B0)]
    B_current <- Assessment@B[length(Assessment@B)]
    if(has_UMSY) {
      SE_UMSY <- ifelse(is.na(Assessment@SE_UMSY[1]), 1e-8, Assessment@SE_UMSY)
      UMSY_vector <- trlnorm(reps, Assessment@UMSY, SE_UMSY)
      TAC <- UMSY_vector * B_current
    }
    if(has_FMSY) {
      SE_FMSY <- ifelse(is.na(Assessment@SE_FMSY[1]), 1e-8, Assessment@SE_FMSY)
      FMSY_vector <- trlnorm(reps, Assessment@FMSY, SE_FMSY)
      TAC <- FMSY_vector * B_current
    }
    ramp <- HCRlin(B_B0, 0.1, 0.4)
    TAC <- ramp * TAC
  } else {
    TAC <- rep(NA, reps) # Missing estimates for HCR. Previous TAC recommendation will be used.
  }

  Rec <- new("Rec")
  Rec@TAC <- TACfilter(TAC)
  return(Rec)
}
class(HCR40_10) <- "HCR"



#' A 60-20 harvest control rule.
#'
#' A output control rule with a ramp that reduces the
#' TAC recommendation (FMSY * Biomass) when biomass depletion (B/B0) is less than 60%.
#' The TAC reduction is linearly reduced to 0 when depletion is less than 20%.
#'
#' @param Assessment An object of class \linkS4class{Assessment}.
#' @param reps The number of stochastic samples of the TAC recommendation.
#' @param ... Miscellaneous arguments (not currently used).
#' @return An object of class \linkS4class{Rec} with the TAC recommendation.
#' @author Q. Huynh
#' @references
#' Punt, A. E, Dorn, M. W., and Haltuch, M. A. 2008. Evaluation of threshold management strategies
#' for groundfish off the U.S. West Coast. Fisheries Research 94:251-266.
#' @seealso \link{HCR40_10} \link{HCRlin} \link{make_MP}
#' @examples
#' # 60-20 linear ramp
#' Brel <- seq(0, 1, length.out = 200)
#' plot(Brel, HCRlin(Brel, 0.2, 0.6), xlab = "Estimated B/B0", ylab = "TAC adjustement factor",
#' main = "A 40-10 harvest control rule", type = 'l', col = 'blue')
#' abline(v = c(0.1,0.4), col = "red", lty = 2)
#'
#' # create an MP to run in closed-loop MSE
#' DD_60_20 <- make_MP(DD_TMB, HCR60_20)
#' \dontrun{
#' myOM <- DLMtool::runMSE(DLMtool::testOM, MPs = c("FMSYref", "DD_60_20"))
#' }
#' @export
HCR60_20 <- function(Assessment, reps = 1, ...) {
  has_UMSY <- has_FMSY <- has_B <- has_B_B0 <- FALSE
  if(!is.na(Assessment@UMSY[1])) has_UMSY <- TRUE
  if(!is.na(Assessment@FMSY[1])) has_FMSY <- TRUE
  if(length(Assessment@B) > 0) has_B <- TRUE
  if(length(Assessment@B_B0) > 0) has_B_B0 <- TRUE

  if((has_UMSY || has_FMSY) && has_B && has_B_B0) {
    B_B0 <- Assessment@B_B0[length(Assessment@B_B0)]
    B_current <- Assessment@B[length(Assessment@B)]
    if(has_UMSY) {
      SE_UMSY <- ifelse(is.na(Assessment@SE_UMSY[1]), 1e-8, Assessment@SE_UMSY)
      UMSY_vector <- trlnorm(reps, Assessment@UMSY, SE_UMSY)
      TAC <- UMSY_vector * B_current
    }
    if(has_FMSY) {
      SE_FMSY <- ifelse(is.na(Assessment@SE_FMSY[1]), 1e-8, Assessment@SE_FMSY)
      FMSY_vector <- trlnorm(reps, Assessment@FMSY, SE_FMSY)
      TAC <- FMSY_vector * B_current
    }
    ramp <- HCRlin(B_B0, 0.2, 0.6)
    TAC <- ramp * TAC
  } else {
    TAC <- rep(NA, reps) # Missing estimates for HCR. Previous TAC recommendation will be used.
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


