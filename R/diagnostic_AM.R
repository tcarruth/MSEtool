
#' diagnostic_AM (diagnostic of Assessments in MSE): Did Assess models converge during MSE?
#'
#' Diagnostic check for convergence of Assess models converge during MSE.
#' Assess models write output to the DLMenv environement if their argument
#' \code{diagnostic = TRUE}.
#'
#' @param MSE An object of class MSE created by \code{\link[DLMtool]{runMSE}}.
#' @param DLMenv The name of the environment that contains the Assessment output
#' generated during the MSE.
#' @param gradient_threshold The value of the maximum gradient magnitude below which the
#' model is considered to have converged.
#' @param figure Logical, whether a figure will be drawn.
#' @return A data frame with diagnostic information for the assesssment-based MPs. If \code{figure = TRUE},
#' a set of figures: traffic light (red/green) plots indicating whether model converged,
#' according to \code{convergence} code in list returned by \code{\link[stats]{nlminb}}
#' the Hessian matrix is positive-definite, according to \code{pdHess} in
#' list returned by \code{\link[TMB]{sdreport}}, and the maximum gradient magnitude is
#' below \code{gradient_threshold}. Model run time during all simulations if the
#' model was timed in the management procedure.
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' DD_MSY <- makeMP(DD_TMB, HCR_MSY, diagnostic = TRUE) # is FALSE by default
#' show(DD_MSY)
#' myMSE <- DLMtool::runMSE(DLMtool::testOM, MPs = "DD_MSY")
#' diagnostic_AM(myMSE)
#'
#' ls(DLMtool::DLMenv) # Assessment output and diagnostics are located here
#'
#' # Save to disk for future use. File may be very large due to size of DLMenv!
#' save(myMSE, DLMenv, myMSE_hist, file = "DLMenv.RData")
#' }
#' @importFrom graphics layout
#' @seealso \link{retrospective_AM}
#' @export
diagnostic_AM <- function(MSE, DLMenv = DLMtool::DLMenv, gradient_threshold = 0.1, figure = TRUE) {
  if(length(ls(DLMenv)) == 0) stop("Nothing found in DLMenv.")

  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5, 4, 1, 1), oma = c(0, 0, 8, 0))
    layout(matrix(c(1, 2, 3, 4, 4, 5), ncol = 3, byrow = TRUE))
  }

  env_objects <- ls(DLMenv)
  MPs <- MSE@MPs
  MPs_in_env <- vapply(MPs, function(x) any(grepl(x, env_objects)), logical(1))
  MPs <- MPs[MPs_in_env]

  get_code <- function(x, y) {
    res <- x[names(x) == y]
    do.call(c, res)
  }

  res_df <- matrix(NA, ncol = length(MPs), nrow = 5)

  message(paste("Creating plots for MP:", paste(MPs, collapse = " ")))
  for(i in 1:length(MPs)) {
    objects_vec <- paste0(c("Assessment_report_", "diagnostic_"), MPs[i])
    objects <- mget(objects_vec, envir = DLMenv, ifnotfound = list(NULL))

    Assessment_report <- objects[[1]]
    diagnostic <- objects[[2]]

    if(!is.null(diagnostic)) {
      convergence_code <- lapply(diagnostic, get_code, y = "conv")
      convergence_code <- do.call(rbind, convergence_code)
      hessian_code <- lapply(diagnostic, get_code, y = "hess")
      hessian_code <- do.call(rbind, hessian_code)
      max_gr <- lapply(diagnostic, get_code, y = "maxgrad")
      max_gr <- do.call(rbind, max_gr)
      opt_time <- lapply(diagnostic, get_code, y = "timing")
      opt_time <- do.call(rbind, opt_time)

      res_df[, i] <- c(100 * (1 - sum(convergence_code)/length(convergence_code)),
                       100 * sum(hessian_code)/length(hessian_code),
                       100 * sum(abs(max_gr) <= gradient_threshold)/length(max_gr),
                       median(opt_time), mean(opt_time))

      if(figure) {
        plot_convergence(convergence_code, "converge")
        plot_convergence(hessian_code, "hessian")
        plot_max_gr(max_gr, gradient_threshold)
        plot_time(opt_time, "line")
        plot_time(opt_time, "hist")
        title(paste(MPs[i], "management procedure"), outer = TRUE)
      }
    }
  }

  res_df <- as.data.frame(round(res_df, 2))
  colnames(res_df) <- MPs
  out_df <- data.frame(Description <- c("Percent convergence", "Percent positive-definite Hessian",
                                        paste0("Percent max. gradient <= ", gradient_threshold),
                                        "Median runtime (seconds)", "Mean runtime (seconds)"))

  return(cbind(out_df, res_df))
}

plot_convergence <- function(convergence_code, plot_type = c('converge', 'hessian')) {
  nsim <- nrow(convergence_code)
  npro <- ncol(convergence_code)

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
      if(plot_type == "converge") color <- ifelse(convergence_code[i, j] == 0, 'green', 'red')
      if(plot_type == "hessian") color <- ifelse(convergence_code[i, j], 'green', 'red')
      xcoord <- rep(c(j - 0.5, j + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}


plot_time <- function(timing, plot_type = c('line', 'hist')) {
  plot_type <- match.arg(plot_type)
  time_range <- range(timing)
  if(max(time_range) > 60) {
    units <- "minutes"
    denom <- 60
  } else {
    units <- "seconds"
    denom <- 1
  }

  if(plot_type == "hist") {
    res <- timing/denom
    hist(res, main = "", xlab = paste0("Run time (", units, ")"))
    mtext(paste("Median:", round(median(res), 2), units))
  }

  if(plot_type == "line") {
    nsim <- nrow(timing)
    npro <- ncol(timing)

    color.vec <- gplots::rich.colors(nsim)
    plot(x = NULL, y = NULL, xlim = c(1, npro),
         ylim = c(0.9, 1.1) * time_range/denom,
         xlab = "MSE projection forward in time", ylab = paste0("Time (", units, ")"))
    mtext('Model runtime during MSE')

    for(i in 1:nsim) {
      points(1:npro, timing[i, ]/denom, col = color.vec[i], typ = 'l')
      text(1:npro, timing[i, ]/denom, labels = i, col = color.vec[i])
    }
  }

  return(invisible())
}


plot_max_gr <- function(max_gr, threshold = 1) {
  nsim <- nrow(max_gr)
  npro <- ncol(max_gr)

  plot(x = NULL, y = NULL, xlim = c(0.5, npro+0.5), ylim = c(0.5, nsim+0.5),
       xlab = "MSE projection forward in time", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  mtext(paste0('Maximum gradient of likelihood \n (green: <= ', threshold, ', red: > ', threshold, ')'))

  for(i in 1:nsim) {
    for(j in 1:npro) {
      color <- ifelse(max_gr[i, j] < threshold, 'green', 'red')
      xcoord <- rep(c(j - 0.5, j + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
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
  diagnostic <- get0(paste0("diagnostic_", MP), envir = DLMenv,
                     ifnotfound = vector("list", nsim))

  # Update reporting objects
  if(inherits(Assessment, "Assessment")) {
    conv <- ifelse(is.character(Assessment@opt), 1L, Assessment@opt$convergence)
    hess <- ifelse(is.character(Assessment@SD), FALSE, Assessment@SD$pdHess)
    maxgrad <- ifelse(is.character(Assessment@SD), 1e10, max(abs(Assessment@SD$gradient.fixed)))
    dg <- list(conv = conv, hess = hess, maxgrad = maxgrad)

    # Remove some objects to save memory/disk space
    #Assessment@obj <- Assessment@info <- Assessment@TMB_report <- list()
    if(hess) Assessment@obj <- list()
    Assessment@info <- Assessment@TMB_report <- list()
    Assessment@SD <- ""
    Assessment@Data <- new("Data", stock = "MSE")
  } else {
    dg <- list(conv = 1L, hess = FALSE, maxgrad = 1e10)
  }

  dg$timing <- ifelse(is.null(opt_timing), NA, as.numeric(opt_timing))

  Assessment_report[[x]] <- c(Assessment_report[[x]], Assessment)
  diagnostic[[x]] <- c(diagnostic[[x]], dg)

  # Assign report objects to DLMenv
  assign(paste0("Assessment_report_", MP), Assessment_report, envir = DLMenv)
  assign(paste0("diagnostic_", MP), diagnostic, envir = DLMenv)

  return(invisible())
}




