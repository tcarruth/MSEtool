
summary_cDD_SS <- function(Assessment) summary_cDD(Assessment, TRUE)

rmd_cDD_SS <- function(Assessment) rmd_cDD(Assessment, TRUE)

#' @importFrom reshape2 acast
profile_likelihood_cDD_SS <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  if(!"R0" %in% names(dots) && !"h" %in% names(dots)) stop("Sequence of neither R0 nor h was not found. See help file.")
  if(!is.null(dots$R0)) R0 <- dots$R0 else {
    R0 <- Assessment@R0
    profile_par <- "h"
  }
  if(!is.null(dots$h)) h <- dots$h else {
    h <- Assessment@h
    profile_par <- "R0"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(R0 = R0, h = h)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$log_R0 <- log(profile_grid[i, 1] * Assessment@info$rescale)
    if(Assessment@info$data$SR_type == "BH") {
      params$transformed_h <- logit((profile_grid[i, 2] - 0.2)/0.8)
    } else {
      params$transformed_h <- log(profile_grid[i, 2] - 0.2)
    }

    if(length(Assessment@opt$par) == 1) { # R0 is the only estimated parameter
      if(!joint_profile && profile_par == "R0") {
        nll <- Assessment@obj$fn(params$log_R0)
      } else {
        obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                          DLL = "MSEtool", silent = TRUE)
        opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
        if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
      }
    } else if(length(Assessment@opt$par) == 2 && all(names(Assessment@opt$par) == c("log_R0", "transformed_h"))) {
      if(joint_profile) {
        nll <- Assessment@obj$fn(c(params$log_R0, params$transformed_h))
      } else {
        if(profile_par == "R0") map$log_R0 <- factor(NA) else map$transformed_h <- factor(NA)
        obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                          DLL = "MSEtool", silent = TRUE)
        opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
        if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
      }
    } else {
      map$log_R0 <- map$transformed_h <- factor(NA)
      obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, random = Assessment@obj$env$random,
                        DLL = "MSEtool", silent = TRUE)
      opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]
      if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    }
    if(!exists("nll")) nll <- NA
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

  if(figure) {
    R0.MLE <- Assessment@R0
    h.MLE <- Assessment@h
    if(joint_profile) {
      z.mat <- acast(profile_grid, list("h", "R0"), value.var = "nll")
      contour(x = h, y = R0, z = z.mat, xlab = "Steepness", ylab = expression(R[0]), nlevels = 20)
      points(h.MLE, R0.MLE, col = "red", cex = 1.5, pch = 16)
    } else {
      if(profile_par == "R0") xlab <- expression(R[0]) else xlab <- "Steepness"
      plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
    }

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      if(joint_profile) {
        z.mat <- acast(profile_grid, list("h", "R0"), value.var = "nll")
        contour(x = h, y = R0, z = z.mat, xlab = "Steepness", ylab = expression(R[0]), nlevels = 20)
        points(h.MLE, R0.MLE, col = "red", cex = 1.5, pch = 16)
        msg <- "Joint profile likelihood of R0 and h. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate."
      } else {
        if(profile_par == "R0") xlab <- expression(R[0]) else xlab <- "Steepness"
        plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
        msg <- paste0("Profile likelihood of ", profile_par, ". Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
      }
      dev.off()

      profile.file.caption <- c("profile_likelihood.png", msg)

      html_report(plot.dir, model = "Continuous Delay-Differential (State-Space)", captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))

    }

  }
  return(profile_grid)
}

retrospective_cDD_SS <- function(Assessment, nyr) retrospective_cDD(Assessment, nyr, TRUE)

