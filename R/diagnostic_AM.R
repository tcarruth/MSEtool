
#' Preliminary Assessments in MSE
#'
#' Evaluates the likely performance of Assessment models in the operating model. This function will
#' apply the assessment model for Data generated during the historical period of the MSE, and report
#' the convergence rate for the model and total time elapsed in running the assessments.
#'
#' @param x Either a \code{Hist}, \code{Data} or \code{OM} object.
#' @param Assess An Assess function of class \code{Assess}.
#' @param ncpus Numeric, the number of CPUs to run the Assessment model (will run in parallel if greater than 1).
#' @param ... Arguments to be passed to \code{Assess}, e.g., model configurations.
#'
#' @return Returns invisibly a list of \linkS4class{Assessment} objects of length \code{OM@@nsim}. Messages via console.
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' prelim_AM(DLMtool::testOM, DD_TMB)
#' }
#' @export
prelim_AM <- function(x, Assess, ncpus = 1, ...) {

  # is Hist?
  if(is.list(x) && any(names(x) == "Data")) {
    Data <- x$Data
  } else
    if(inherits(x, "OM")) {
      message("Generating Hist object from OM object via runMSE:")
      runHist <- runMSE(x, Hist = TRUE)
      Data <- runHist$Data
    } else
      if(inherits(x, "Data")) {
        Data <- x
      } else {
        stop("x does not appear to be either a Hist, Data, or OM object.")
      }

  snowfall::sfInit(parallel = ncpus > 1, cpus = ncpus)
  snowfall::sfLibrary("MSEtool", character.only = TRUE)
  nsim <- nrow(Data@Cat)
  message(paste0("Running ", deparse(substitute(Assess)), " with ", nsim, " simulations for ", deparse(substitute(x)), "."))
  dots <- list(...)
  if(length(dots) > 0) message(paste0("\nAdditional arguments to be provided to ", deparse(substitute(Assess)), ":\n", paste(names(dots), collapse = "\n")))
  Assess <- match.fun(Assess)
  if(!inherits(Assess, "Assess")) stop(paste(deparse(substitute(Assess))), "does not appear to be an Assess function.")

  if(snowfall::sfParallel()) snowfall::sfExport(list = c("Assess", "Data"))
  timing <- proc.time()
  res <- snowfall::sfLapply(1:nsim, Assess, Data = Data, ...)
  timing2 <- (proc.time() - timing)[3]
  message("Assessments complete.")

  snowfall::sfStop()

  message(paste0("Total time to run ", nsim, " assessments: ", round(timing2, 1), " seconds"))

  nonconv <- !vapply(res, getElement, logical(1), "conv")
  message(paste0(sum(nonconv), " of ", nsim, " simulations (", round(100 *sum(nonconv)/nsim, 1), "%) failed to converge."))
  if(sum(nonconv > 0)) message(paste("See simulation number:", paste(which(nonconv), collapse = " ")))

  return(invisible(res))
}





#' diagnostic_AM (diagnostic of Assessments in MSE): Did Assess models converge during MSE?
#'
#' Diagnostic check for convergence of Assess models during MSE.
#' Assess models write output to the DLMenv environment if the MP was created with \link{make_MP}
#' with argument \code{diagnostic = TRUE}. This function summarizes and plots the diagnostic information.
#'
#' @param MSE An object of class MSE created by \code{\link[DLMtool]{runMSE}}. If no MSE object
#' is available, use argument \code{MP} instead.
#' @param DLMenv The name of the environment that contains the Assessment output
#' generated during the MSE.
#' @param MP A character vector of MPs with assessment models.
#' @param gradient_threshold The value of the maximum gradient magnitude below which the
#' model is considered to have converged.
#' @param figure Logical, whether a figure will be drawn.
#' @return A matrix with diagnostic information for the assesssment-based MPs. If \code{figure = TRUE},
#' a set of figures: traffic light (red/green) plots indicating whether model converged,
#' according to \code{convergence} code in list returned by \code{\link[stats]{nlminb}}
#' the Hessian matrix is positive-definite, according to \code{pdHess} in
#' list returned by \code{\link[TMB]{sdreport}}, and the maximum gradient magnitude is
#' below \code{gradient_threshold}. Also includes model runtime, number of optimization iterations, and
#' number of function evaluations of the assessment model during each application of the management procedure.
#' @author Q. Huynh
#' @examples
#' \dontrun{
#' DD_MSY <- makeMP(DD_TMB, HCR_MSY, diagnostic = "min")
#' show(DD_MSY)
#' myMSE <- runMSE(DLMtool::testOM, MPs = "DD_MSY")
#' diagnostic_AM(myMSE)
#'
#' # If MSE object is not available (e.g. runMSE crashed), use MP argument instead.
#' diagnostic_AM(MP = "DD_MSY")
#'
#' ls(DLMtool::DLMenv) # Assessment output and diagnostics are located here
#'
#' # Save to disk for future use. File may be very large due to size of DLMenv!
#' save(myMSE, DLMenv, file = "DLMenv.RData")
#' }
#' @importFrom graphics layout
#' @seealso \link{retrospective_AM}
#' @export
diagnostic_AM <- function(MSE = NULL, DLMenv = DLMtool::DLMenv, MP = NULL, gradient_threshold = 0.1, figure = TRUE) {
  if(length(ls(DLMenv)) == 0) stop("Nothing found in DLMenv.")

  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(layout(matrix(1)))
    on.exit(par(old_par), add = TRUE)

    par(mar = c(5, 4, 1, 1), oma = c(0, 0, 8, 0))
  }

  env_objects <- ls(DLMenv)
  if(!is.null(MSE)) {
    MPs <- MSE@MPs
  } else if(!is.null(MP)) {
    MPs <- MP
  } else stop("No MSE object or character vector of MP was provided.")
  MPs_in_env <- vapply(MPs, function(x) any(grepl(x, env_objects)), logical(1))
  MPs <- MPs[MPs_in_env]

  get_code <- function(x, y) vapply(x, getElement, numeric(1), y)

  res_mat <- matrix(NA, ncol = length(MPs), nrow = 5)

  message(paste0("Creating plots for MP:\n", paste(MPs, collapse = "\n")))
  for(i in 1:length(MPs)) {
    objects_vec <- paste0(c("Assessment_report_", "diagnostic_"), MPs[i])
    objects <- mget(objects_vec, envir = DLMenv, ifnotfound = list(NULL))

    Assessment_report <- objects[[1]]
    diagnostic <- objects[[2]]

    if(!is.null(diagnostic)) {
      convergence_code <- lapply(diagnostic, get_code, y = "conv")
      hessian_code <- lapply(diagnostic, get_code, y = "hess")
      max_gr <- lapply(diagnostic, get_code, y = "maxgrad")
      iter <- lapply(diagnostic, get_code, y = "iter")
      fn_eval <- lapply(diagnostic, get_code, y = "fn_eval")
      Year <- lapply(diagnostic, get_code, y = "Year")

      if(figure) {
        layout(matrix(c(1, 2, 3, 4, 4, 5), ncol = 3, byrow = TRUE))

        plot_convergence(convergence_code, Year, "converge")
        plot_convergence(hessian_code, Year, "hessian")
        plot_max_gr(max_gr, Year, gradient_threshold)

        plot_iter(iter, Year, "line", "Optimization iterations")
        plot_iter(fn_eval, Year, "hist", "Function evaluations")

        title(paste(MPs[i], "management procedure"), outer = TRUE)
      }

      convergence_code <- do.call(c, convergence_code)
      hessian_code <- do.call(c, hessian_code)
      max_gr <- do.call(c, max_gr)
      iter <- do.call(c, iter)
      fn_eval <- do.call(c, fn_eval)

      res_mat[, i] <- c(100 * (1 - sum(convergence_code)/length(convergence_code)),
                        100 * sum(hessian_code)/length(hessian_code),
                        100 * sum(abs(max_gr) <= gradient_threshold)/length(max_gr),
                        median(iter, na.rm = TRUE), median(fn_eval, na.rm = TRUE))
    }
  }

  res_mat <- round(res_mat, 2)
  colnames(res_mat) <- MPs
  rownames(res_mat) <- c("Percent convergence", "Percent positive-definite Hessian",
                         paste0("Percent max. gradient <= ", gradient_threshold),
                         "Median iterations", "Median function evaluations")
  return(res_mat)
}

plot_convergence <- function(convergence_code, Year, plot_type = c('converge', 'hessian')) {
  nsim <- length(convergence_code)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE year", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  if(plot_type == "converge") {
    mtext('Did nlminb optimizer \n report model convergence? \n (green: yes, red: no)')
  }
  if(plot_type == "hessian") {
    mtext('Is Hessian matrix \n positive definite? \n (green: yes, red: no)')
  }
  for(i in 1:nsim) {
    for(j in 1:length(convergence_code[[i]])) {
      if(plot_type == "converge") color <- ifelse(convergence_code[[i]][j] == 0, 'green', 'red')
      if(plot_type == "hessian") color <- ifelse(convergence_code[[i]][j], 'green', 'red')
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}


plot_iter <- function(x, Year, plot_type = c('line', 'hist'), lab = c("Optimization iterations", "Function evaluations")) {
  plot_type <- match.arg(plot_type)
  lab <- match.arg(lab)

  if(plot_type == "hist") {
    hist(do.call(c, x), main = "", xlab = lab)
    mtext(paste("Median:", round(median(do.call(c, x)), 2)))
  }

  if(plot_type == "line") {
    nsim <- length(x)
    yr_range <- range(do.call(c, Year))

    color.vec <- gplots::rich.colors(nsim)
    plot(x = NULL, y = NULL, xlim = yr_range,
         ylim = c(0.9, 1.1) * range(do.call(c, x), na.rm = TRUE),
         xlab = "MSE Year", ylab = lab)
    mtext(lab)

    for(i in 1:nsim) {
      points(Year[[i]], x[[i]], col = color.vec[i], typ = 'l')
      text(Year[[i]], x[[i]], labels = i, col = color.vec[i])
    }
  }

  return(invisible())
}

plot_max_gr <- function(max_gr, Year, threshold = 1) {
  nsim <- length(max_gr)
  yr_range <- range(do.call(c, Year))

  plot(x = NULL, y = NULL, xlim = c(-0.5, 0.5) + yr_range, ylim = c(0.5, nsim + 0.5),
       xlab = "MSE Year", ylab = "Simulation #", xaxs = "i", yaxs = "i")
  mtext(paste0('Maximum gradient of likelihood \n (green: <= ', threshold, ', red: > ', threshold, ')'))

  for(i in 1:nsim) {
    for(j in 1:length(max_gr[[i]])) {
      color <- ifelse(max_gr[[i]][j] < threshold, 'green', 'red')
      xcoord <- rep(c(Year[[i]][j] - 0.5, Year[[i]][j] + 0.5), each = 2)
      ycoord <- c(i - 0.5, i + 0.5, i + 0.5, i - 0.5)
      polygon(x = xcoord, y = ycoord, border = 'black', col = color)
    }
  }
  return(invisible())
}



# Prints Assessment during MSE output to DLMtool::DLMenv
# Call in MP created by make_MP when diagnostic = "min" or "full"
Assess_diagnostic <- function(DLMenv = DLMtool::DLMenv, include_assessment = TRUE) {
  # Get Assessment objects
  report_obj <- mget(c("x", "do_Assessment", "Data"), envir = parent.frame(),
                     ifnotfound = list(NULL))
  x <- report_obj$x
  Assessment <- report_obj$do_Assessment
  Data <- report_obj$Data

  # Check if this function is called within an MSE
  runMSE.called <- any(vapply(sys.calls(), function(x) any(grepl("runMSE", x)), logical(1)))
  if(!runMSE.called) return(invisible()) # exit, no call to runMSE found in sys.calls()

  # Assume passed checks
  MP <- eval.parent(expression(MPs[mp]), n = 3)
  nsim <- nrow(Data@OM)
  Year <- Data@Year[length(Data@Year)]

  # Update reporting objects
  if(inherits(Assessment, "Assessment")) {
    conv <- ifelse(is.character(Assessment@opt), 1L, Assessment@opt$convergence)
    hess <- ifelse(is.character(Assessment@SD), FALSE, Assessment@SD$pdHess)
    maxgrad <- ifelse(is.character(Assessment@SD), 1e10, max(abs(Assessment@SD$gradient.fixed)))
    iter <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$iterations)
    fn_eval <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$evaluations[1])

    dg <- list(conv = conv, hess = hess, maxgrad = maxgrad, iter = iter, fn_eval = fn_eval, Year = Year)

  } else {
    dg <- list(conv = 1L, hess = FALSE, maxgrad = 1e10, iter = NA, fn_eval = NA, Year = Year)
  }

  # Assign report objects to DLMenv
  diagnostic <- get0(paste0("diagnostic_", MP), envir = DLMenv,
                     ifnotfound = vector("list", nsim))
  len_diag <- length(diagnostic[[x]])
  diagnostic[[x]][[len_diag + 1]] <- dg
  assign(paste0("diagnostic_", MP), diagnostic, envir = DLMenv)

  if(include_assessment) {
    # Remove some objects to save memory/disk space
    if(inherits(Assessment, "Assessment")) {
      if(dg$hess) Assessment@obj <- list()
      Assessment@info <- Assessment@TMB_report <- list()
      Assessment@SD <- ""
    }

    Assessment_report <- get0(paste0("Assessment_report_", MP), envir = DLMenv,
                              ifnotfound = vector("list", nsim))
    len_Assess <- length(Assessment_report[[x]])
    if(len_Assess > 0) {
      Assessment_report[[x]][[len_Assess + 1]] <- Assessment
    } else Assessment_report[[x]] <- c(Assessment_report[[x]], Assessment)
    assign(paste0("Assessment_report_", MP), Assessment_report, envir = DLMenv)
  }

  return(invisible())
}


# Assign report objects to DLMenv
#if(!exists(paste0("diagnostic_", MP), envir = DLMenv)) {
#  assign(paste0("diagnostic_", MP), vector("list", nsim), envir = DLMenv)
#}
#dg_exp <- bquote({
#  len_diag <- length(.(as.symbol(paste0("diagnostic_", MP)))[[.(x)]])
#  .(as.symbol(paste0("diagnostic_", MP)))[[.(x)]][[len_diag + 1]] <- get("dg", envir = .(sys.frames()[[length(sys.frames())]]))
#})
#eval(dg_exp, envir = DLMenv)

#if(!exists(paste0("Assessment_report_", MP), envir = DLMenv)) {
#  assign(paste0("Assessment_report_", MP), vector("list", nsim), envir = DLMenv)
#}
#Assess_exp <- bquote({
#  len_Assess <- length(.(as.symbol(paste0("Assessment_report_", MP)))[[.(x)]])
#  if(len_Assess > 0) {
#    .(as.symbol(paste0("Assessment_report_", MP)))[[.(x)]][[len_diag + 1]] <- get("Assessment", envir = .(sys.frames()[[length(sys.frames())]]))
#  } else {
#    .(as.symbol(paste0("Assessment_report_", MP)))[[.(x)]] <- c(.(as.symbol(paste0("Assessment_report_", MP)))[[.(x)]],
#                                                                get("Assessment", envir = .(sys.frames()[[length(sys.frames())]])))
#  }
#})
#eval(Assess_exp, envir = DLMenv)


