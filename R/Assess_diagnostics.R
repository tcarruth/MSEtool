
#' Do Assess models converge during MSE?
#'
#' Diagnostic check for convergence of Assess models converge during MSE.
#' Assess models write output to the DLMenv environement if their argument
#' \code{diagnostic = TRUE}.
#' @param MSE An object of class MSE created by \code{\link[DLMtool]{runMSE}}.
#' @param DLMenv The name of the environment that contains the Assessment output
#' generated during the MSE.
#' @param gradient_threshold The value of the maximum gradient magnitude below which the
#' model is considered to have converged.
#' @return Traffic light (red/green) plots indicating whether model converged,
#' according to \code{convergence} code in list returned by \code{\link[stats]{nlminb}}
#' the Hessian matrix is positive-definite, according to \code{pdHess} in
#' list returned by \code{\link[TMB]{sdreport}}, and the maximum gradient magnitude is
#' below \code{gradient_threshold}. Also, model run time during all simulations if the
#' model was timed in the management procedure.
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' DD_MSY <- makeMP(DD_TMB, HCR_MSY, diagnostic = TRUE) # is TRUE by default
#' show(DD_MSY)
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = "DD_MSY")
#' plot_diagnostics(myMSE)
#'
#' ls(DLMtool::DLMenv) # Assessment output and diagnostics are located here
#' save(myMSE, DLMenv, file = "DLMenv.RData") # Save to disk. Note: very large object!
#' }
#' @importFrom graphics layout
#' @export
plot_diagnostics <- function(MSE, DLMenv = DLMtool::DLMenv, gradient_threshold = 0.1) {
  if(length(ls(DLMenv)) == 0) stop(paste0("Nothing found in ", DLMenv))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  env_objects <- ls(DLMenv)
  MPs <- MSE@MPs
  MPs_in_env <- vapply(MPs, function(x) any(grepl(x, env_objects)), logical(1))
  MPs <- MPs[MPs_in_env]

  get_convergence_code <- function(x) {
    vapply(x, function(y) ifelse(is.character(y@opt), 1L, y@opt$convergence), numeric(1))
  }
  get_hessian_code <- function(x) {
    vapply(x, function(y) ifelse(is.character(y@SD), FALSE, y@SD$pdHess), logical(1))
  }

  get_max_gr <- function(x) {
    vapply(x, function(y) ifelse(is.character(y@SD), 1e10, max(abs(y@SD$gradient.fixed))), numeric(1))
  }

  message(paste("Creating plots for MP:", paste(MPs, collapse = " ")))
  for(i in 1:length(MPs)) {
    objects_vec <- paste0(c("Assessment_report_", "opt_time_"), MPs[i])
    objects <- mget(objects_vec, envir = DLMenv, ifnotfound = list(NULL))

    Assessment_report <- objects[[1]]
    opt_time <- objects[[2]]

    if(!is.null(Assessment_report)) {
      if(!is.null(opt_time)) {
        layout(matrix(c(1, 2, 3, 4, 4, 5), ncol = 3, byrow = TRUE))
      } else par(mfrow = c(1, 3))
      par(oma = c(0, 0, 2, 0))

      convergence_code <- lapply(Assessment_report, get_convergence_code)
      plot_convergence(convergence_code, "converge")

      hessian_code <- lapply(Assessment_report, get_hessian_code)
      plot_convergence(hessian_code, "hessian")

      max_gr <- lapply(Assessment_report, get_max_gr)
      plot_max_gr(max_gr, gradient_threshold)
    }

    if(!is.null(opt_time)) {
      plot_time(opt_time, "line")
      plot_time(opt_time, "hist")
    }
    title(paste(MPs[i], "management procedure"), outer = TRUE)
  }

  return(invisible())
}

plot_convergence <- function(convergence_code, plot_type = c('converge', 'hessian')) {
  nsim <- length(convergence_code)
  npro <- length(convergence_code[[1]])

  plot(x = NULL, y = NULL, xlim = c(0.5, npro+0.5), ylim = c(0.5, nsim+0.5),
       xlab = "MSE projection forward in time", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  if(plot_type == "converge") {
    mtext('Did nlminb optimizer \n report model convergence? \n (green: yes, red: no)')
  }
  if(plot_type == "hessian") {
    mtext('Is Hessian matrix \n positive definite? \n (green: yes, red: no)')
  }
  for(i in 1:nsim) {
    for(j in 1:npro) {
      if(plot_type == "converge") color <- ifelse(convergence_code[[i]][j] == 0, 'green', 'red')
      if(plot_type == "hessian") color <- ifelse(convergence_code[[i]][j], 'green', 'red')
      xcoord <- rep(c(j - 0.5, j + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}


plot_time <- function(timing, plot_type = c('line', 'hist')) {
  plot_type <- match.arg(plot_type)
  time_range <- range(do.call(c, timing))
  if(max(time_range) > 60) {
    units <- "minutes"
    denom <- 60
  } else {
    units <- "seconds"
    denom <- 1
  }

  if(plot_type == "hist") {
    res <- do.call(c, timing)/denom
    hist(res, main = "", xlab = paste0("Run time (", units, ")"))
    mtext(paste("Median:", round(median(res), 2), units))
  }

  if(plot_type == "line") {
    nsim <- length(timing)
    npro <- length(timing[[1]])

    color.vec <- gplots::rich.colors(nsim)
    plot(x = NULL, y = NULL, xlim = c(1, npro),
         ylim = c(0.9, 1.1) * time_range/denom,
         xlab = "MSE projection forward in time", ylab = paste0("Time (", units, ")"))
    mtext('Model runtime during MSE')

    for(i in 1:nsim) {
      points(1:npro, timing[[i]]/denom, col = color.vec[i], typ = 'l')
      text(1:npro, timing[[i]]/denom, labels = i, col = color.vec[i])
      #text(npro, timing[[i]][npro]/denom, labels = i, pos = 3, col = color.vec[i])
    }
  }

  return(invisible())
}


plot_max_gr <- function(max_gr, threshold = 1) {
  nsim <- length(max_gr)
  npro <- length(max_gr[[1]])

  plot(x = NULL, y = NULL, xlim = c(0.5, npro+0.5), ylim = c(0.5, nsim+0.5),
       xlab = "MSE projection forward in time", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  mtext(paste0('Maximum gradient of likelihood \n (green: <= ', threshold, ', red: > ', threshold, ')'))

  for(i in 1:nsim) {
    for(j in 1:npro) {
      color <- ifelse(max_gr[[i]][j] < threshold, 'green', 'red')
      xcoord <- rep(c(j - 0.5, j + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}


retrospective_MSE <- function(MSE, DLMenv = DLMtool::DLMenv, sim = 1, MP,
                              plot_type = c("SSB", "F", "SSB_SSBMSY", "F_FMSY", "SSB_SSB0", "VB")) {
  if(length(ls(DLMenv)) == 0) stop(paste0("Nothing found in ", DLMenv))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfcol = c(2, 3), mar = c(5, 4, 1, 1), oma = c(0, 0, 2.5, 0))

  env_objects <- ls(DLMenv)
  MPs <- MSE@MPs
  MPs_in_env <- vapply(MPs, function(x) any(grepl(x, env_objects)), logical(1))
  MPs <- MPs[MPs_in_env]

  plot_type <- match.arg(plot_type, several.ok = TRUE)
  if(!MP %in% MSE@MPs || !MP %in% MPs) stop(paste(MP, "was not found in MSE object (MSE@MPs)."))
  if(!MP %in% MPs) stop(paste(MP, "was not found in DLMenv."))
  if(length(sim) > 1 || sim > MSE@nsim) stop(paste0(sim, " should be a number between 1 and ", MSE@nsim, "."))
  MPind <- MSE@MPs == MP
  if(sum(MPind) > 1) stop(paste("More than one match to", MP, "was found in MSE object."))

  Assessment_report <- get(paste0("Assessment_report_", MP), envir = DLMenv)[[sim]]
  isSP <- grepl("SP", Assessment_report[[1]]@Model)

  color.vec <- gplots::rich.colors(length(Assessment_report))
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
      #Assess <- lapply(Assessment_report, slot, "F_FMSY")
      AssessU <- lapply(Assessment_report, slot, "U")
      AssessUMSY <- lapply(Assessment_report, slot, "UMSY")
      Assess <- Map(function(x, y) log(1-x)/log(1-y), x = AssessU, y = AssessUMSY)
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
      #Assess <- lapply(Assessment_report, slot, "F")
      AssessU <- lapply(Assessment_report, slot, "U")
      Assess <- lapply(AssessU, function(x) -log(1 - x))
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
      Hist <- apply(MSE@VB_hist, c(1, 3), sum)[sim, ]
      Proj <- MSE@VB[sim, MPind, ]
      Assess <- lapply(Assessment_report, slot, "VB")
    }

    if(plot_type[i] != "MSY") {
      xlimits <- c(1, MSE@nyears + MSE@proyears)
      ylimits <- c(0, 1.1 * max(c(Hist, Proj, do.call(c, Assess))))
      plot(c(Hist, Proj), xlab = "MSE year", ylab = plot_type[i], xlim = xlimits, ylim = ylimits, lwd = 2, typ = 'l')
      for(j in length(Assess):1) {
        if(length(Assess[[j]]) > 0) lines(1:length(Assess[[j]]), Assess[[j]], col = color.vec[j])
      }
      if(plot_type[i] == "SSB_SSBMSY" || plot_type[i] == "F_FMSY") abline(h = 1)
      abline(h = 0, col = "grey")
      abline(v = MSE@nyears, lty = 2)
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


# Prints Assessment during MSE output to DLMtool::DLMenv
# Call in MP created by make_MP when diagnostic = TRUE
Assess_diagnostic <- function(DLMenv = DLMtool::DLMenv) {
  # Get Assessment objects
  report_obj <- mget(c("x", "opt_timing", "do_Assessment"), envir = parent.frame(),
                     ifnotfound = list(NULL))
  x <- report_obj$x
  opt_timing <- report_obj$opt_timing
  Assessment <- report_obj$do_Assessment

  # Check if this function is called within an MSE
  Data <- Assessment@Data
  if(ncol(Data@OM) == 0) return(invisible()) # exit, not an MSE
  runMSE.called <- any(vapply(sys.calls(), function(x) any(grepl("runMSE", x)), logical(1)))
  if(!runMSE.called) return(invisible()) # exit, no call to runMSE found in sys.calls()

  # Assume passed checks
  MP <- eval.parent(expression(MPs[mp]), n = 3)
  nsim <- nrow(Data@OM)

  Assessment_report <- get0(paste0("Assessment_report_", MP), envir = DLMenv,
                            ifnotfound = vector("list", nsim))
  opt_time <- get0(paste0("opt_time_", MP), envir = DLMenv,
                   ifnotfound = vector("list", nsim))

  # Update reporting objects
  if(inherits(Assessment, "Assessment")) { # Remove some objects to save memory/disk space
    Assessment@obj <- Assessment@info <- Assessment@TMB_report <- list()
    Assessment@SD <- ""
    Assessment@Data <- new("Data", stock = "MSE")
  }
  Assessment_report[[x]] <- c(Assessment_report[[x]], Assessment)
  opt_time[[x]] <- c(opt_time[[x]], opt_timing)

  # Assign report objects to DLMenv
  assign(paste0("Assessment_report_", MP), Assessment_report, envir = DLMenv)
  assign(paste0("opt_time_", MP), opt_time, envir = DLMenv)

  return(invisible())
}




