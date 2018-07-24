#' retrospective_AM (retrospective of Assessment model in MSE)
#'
#' Plots the true retrospective of an assessment model during the MSE. A series of time series estimates of SSB, F, and VB
#' are plotted over the course of the MSE are plotted against the operating model (true) values (in black).
#'
#' @param MSE An object of class MSE created by \code{\link[DLMtool]{runMSE}}.
#' @param DLMenv The name of the environment that contains the Assessment output
#' generated during the MSE.
#' @param sim Integer between 1 and MSE@@nsim. The simulation number for which the retrospectives will be plotted.
#' @param MP Character. The name of the management procedure created by \code{\link{make_MP}} containing the asssessment model.
#' @param MSE_Hist Optional. The list containing historical data for the MSE, created by \code{\link[DLMtool]{runMSE}} with argument \code{Hist = TRUE}.
#' Currently only used to plot operating model vulnerable biomass in historical period.
#' @param plot_legend Logical. Whether to plot legend to reference year of assessment in the MSE.
#' @author Q. Huynh
#' @details For assessment models that utilize annual harvest rates (u), the instantaneous fishing mortality rates
#' are obtained as F = -log(1 - u).
#' @note This function only plots retrospectives from a single simulation in the MSE. Results from one figure
#' may not be indicative of general assessment behavior and performance overall.
#'
#' For \link{SP} and \link{SP_SS} assessment models don't model SSB. Instead, the estimated vulnerable biomass is plotted.
#' @return A series of figures for spawning stock biomass
#' (SSB, including absolute magnitude and relative to MSY and virgin), fishing mortality (F, including absolute
#' magnitude and relative to MSY), and vulnerable biomass (VB) estimates over the course of the MSE are plotted
#' against the operating model (true) values (in black).
#' @examples
#' \dontrun{
#' DD_MSY <- makeMP(DD_TMB, HCR_MSY, diagnostic = "full")
#' myMSE_hist <- DLMtool::runMSE(DLMtool::testOM, MPs = "DD_MSY", MSE_Hist = TRUE)
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = "DD_MSY")
#' retrospective_AM(myMSE, sim = 1, MP = "DD_MSY")
#' retrospective_AM(myMSE, sim = 1, MP = "DD_MSY", Hist = myMSE_hist)
#'
#' ls(DLMtool::DLMenv) # Assessment output and diagnostics are located here
#'
#' # Save to disk for future use. File may be very large due to size of DLMenv!
#' save(myMSE, DLMenv, myMSE_hist, file = "DLMenv.RData")
#' }
#'
#' @seealso \link{diagnostic_AM}
#' @export
retrospective_AM <- function(MSE, DLMenv = DLMtool::DLMenv, sim = 1, MP, MSE_Hist = NULL, plot_legend = FALSE) {

  plot_type = c("SSB", "F", "SSB_SSBMSY", "F_FMSY", "SSB_SSB0", "VB")

  if(length(ls(DLMenv)) == 0) stop(paste0("Nothing found in ", DLMenv))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfcol = c(2, 3), mar = c(5, 4, 1, 1), oma = c(0, 0, 2.5, 0))

  env_objects <- ls(DLMenv)
  MPs <- MSE@MPs
  MPs_in_env <- vapply(MPs, function(x) any(grepl(x, env_objects)), logical(1))
  MPs <- MPs[MPs_in_env]

  if(!MP %in% MSE@MPs || !MP %in% MPs) stop(paste(MP, "was not found in MSE object (MSE@MPs)."))
  if(!MP %in% MPs) stop(paste(MP, "was not found in DLMenv."))
  if(length(sim) > 1 || sim > MSE@nsim) stop(paste0(sim, " should be a number between 1 and ", MSE@nsim, "."))
  MPind <- MSE@MPs == MP
  if(sum(MPind) > 1) stop(paste("More than one match to", MP, "was found in MSE object."))

  Assessment_report <- get(paste0("Assessment_report_", MP), envir = DLMenv)[[sim]]
  isSP <- grepl("SP", Assessment_report[[1]]@Model)

  color.vec <- gplots::rich.colors(length(Assessment_report))

  Yr_MSE <- 1:(MSE@nyears + MSE@proyears)

  Assess <- lapply(Assessment_report, slot, "VB")
  End_Assess_Yr <- vapply(Assess, function(x) as.numeric(names(x))[length(x)-1], numeric(1))
  for(i in 1:length(plot_type)) {
    if(plot_type[i] == "SSB_SSBMSY") {
      Hist <- apply(MSE@SSB_hist, c(1, 3), sum)[sim, ]/MSE@OM$SSBMSY[sim]
      Proj <- MSE@B_BMSY[sim, MPind, ]
      if(!isSP) {
        Assess <- lapply(Assessment_report, slot, "SSB_SSBMSY")
      } else {
        Assess <- lapply(Assessment_report, slot, "VB_VBMSY")
      }
    }
    if(plot_type[i] == "F_FMSY") {
      Hist <- apply(MSE@FM_hist, c(1, 3), max)[sim, ]/MSE@OM$FMSY[sim]
      Proj <- MSE@F_FMSY[sim, MPind, ]
      Assess <- lapply(Assessment_report, slot, "F_FMSY")
      if(length(do.call(c, Assess)) == 0) {
        AssessU <- lapply(Assessment_report, slot, "U")
        AssessUMSY <- lapply(Assessment_report, slot, "UMSY")
        Assess <- Map(function(x, y) log(1-x)/log(1-y), x = AssessU, y = AssessUMSY)
      }
    }
    if(plot_type[i] == "SSB") {
      Hist <- apply(MSE@SSB_hist, c(1, 3), sum)[sim, ]
      Proj <- MSE@SSB[sim, MPind, ]
      if(!isSP) {
        Assess <- lapply(Assessment_report, slot, "SSB")
      } else {
        Assess <- lapply(Assessment_report, slot, "VB")
      }
    }
    if(plot_type[i] == "F") {
      Hist <- apply(MSE@FM_hist, c(1, 3), max)[sim, ]
      Proj <- MSE@FM[sim, MPind, ]
      Assess <- lapply(Assessment_report, slot, "FMort")
      if(length(do.call(c, Assess)) == 0) {
        AssessU <- lapply(Assessment_report, slot, "U")
        Assess <- lapply(AssessU, function(x) -log(1 - x))
      }
    }
    if(plot_type[i] == "SSB_SSB0") {
      Hist <- apply(MSE@SSB_hist, c(1, 3), sum)[sim, ]/MSE@OM$SSB0[sim]
      Proj <- MSE@SSB[sim, MPind, ]/MSE@OM$SSB0[sim]
      if(!isSP) {
        Assess <- lapply(Assessment_report, slot, "SSB_SSB0")
      } else {
        Assess <- lapply(Assessment_report, slot, "VB_VB0")
      }
    }
    if(plot_type[i] == "VB") {
      if(!is.null(MSE_Hist)) {
        Hist <- MSE_Hist$TSdata$VB[, sim]
      } else {
        Hist <- rep(NA, MSE@nyears)
        message("Provide Hist object in order to plot simulated vulnerable biomass.")
      }
      Proj <- MSE@VB[sim, MPind, ]
      Assess <- lapply(Assessment_report, slot, "VB")
    }

    if(plot_type[i] == "Recruit") {
      if(!is.null(MSE_Hist)) {
        Hist <- MSE_Hist$AtAge$Nage[sim, 1, ]
      } else {
        Hist <- rep(NA, MSE@nyears)
        message("Provide Hist object in order to plot simulated recruitment.")
      }
      Proj <- rep(NA, MSE@proyears)
      Assess <- lapply(Assessment_report, slot, "R")
    }

    if(plot_type[i] != "MSY") {
      Assess_Yr <- lapply(Assess, function(x) as.numeric(names(x)))
      xlimits <- c(1, MSE@nyears + MSE@proyears)

      converged_assessments <- vapply(Assessment_report, slot, logical(1), "conv")
      ylimits <- c(0, 1.1 * max(c(Hist, Proj, do.call(c, Assess[converged_assessments])), na.rm = TRUE))

      if(all(!is.na(ylimits))) {
        plot(Yr_MSE, c(Hist, Proj), xlab = "MSE year", ylab = plot_type[i], xlim = xlimits, ylim = ylimits, lwd = 2, typ = 'l')
        for(j in length(Assess):1) {
          if(Assessment_report[[j]]@conv && length(Assess[[j]]) > 0) lines(Assess_Yr[[j]], Assess[[j]], col = color.vec[j])
        }
        if(plot_type[i] == "SSB_SSBMSY" || plot_type[i] == "F_FMSY") abline(h = 1)
        abline(h = 0, col = "grey")
        abline(v = MSE@nyears, lty = 2)
        if(plot_legend && i == 1) legend("topleft", c("OM", End_Assess_Yr), col = c("black", color.vec), lwd = c(2, rep(1, length(Assessment_report))))
      } else {
        message(paste0("Skipped plot for ", plot_type[i], "."))
      }
    }

    if(plot_type[i] == "MSY") {
      Hist <- c(MSE@OM$FMSY[sim], MSE@OM$MSY[sim])
      Assess <- vapply(Assessment_report, function(x) c(-log(1 - slot(x, "UMSY")), slot(x, "MSY")), numeric(2))

      xc <- c(Hist[1], Assess[1, ])
      xc <- xc[xc < 2]
      xlimits <- c(max(0, mean(xc) - 2*sd(xc)), mean(xc + 2*sd(xc))) # FMSY
      yc <- c(Hist[2], Assess[2, ])
      ylimits <- c(max(0, mean(yc) - 2*sd(yc)), mean(yc + 2*sd(yc))) # MSY
      plot(x = Assess[1, ], y = Assess[2, ], xlim = xlimits, ylim = ylimits,
           xlab = "FMSY estimates", ylab = "MSY estimates", col = color.vec, pch = 16,
           cex = 2)
      j <- ncol(Assess) - 1
      arrows(x0 = Assess[1, 1:j], y0 = Assess[2, 1:j], x1 = Assess[1, 2:(j+1)], y1 = Assess[2, 2:(j+1)],
             length = 0.1)
      points(Hist[1], Hist[2], col = "grey", cex = 2, pch = 0)
    }

  }
  title(paste0(MP, " management procedure \n Simulation #", sim), outer = TRUE)
  invisible()
}
