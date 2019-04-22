summary_SP_SS <- function(Assessment) summary_SP(Assessment, TRUE)

rmd_SP_SS <- function(Assessment) rmd_SP(Assessment, TRUE)




#' @importFrom reshape2 acast
profile_likelihood_SP_SS <- function(Assessment, figure = TRUE, save_figure = TRUE, save_dir = tempdir(), ...) {
  dots <- list(...)
  dots <- list(...)
  if(!"FMSY" %in% names(dots) && !"MSY" %in% names(dots)) stop("Sequence of neither FMSY nor MSY was found. See help file.")
  if(!is.null(dots$FMSY)) FMSY <- dots$FMSY else {
    FMSY <- Assessment@FMSY
    profile_par <- "MSY"
  }
  if(!is.null(dots$MSY)) MSY <- dots$MSY else {
    MSY <- Assessment@MSY
    profile_par <- "FMSY"
  }

  map <- Assessment@obj$env$map
  params <- Assessment@info$params

  profile_grid <- expand.grid(FMSY = FMSY, MSY = MSY)
  joint_profile <- !exists("profile_par")

  profile_fn <- function(i, Assessment, params, map) {
    params$log_FMSY <- log(profile_grid[i, 1])
    params$log_MSY <- log(profile_grid[i, 2] * Assessment@info$rescale)

    if(joint_profile) {
	  map$log_MSY <- map$log_FMSY <- factor(NA)
	} else {
	  if(profile_par == "MSY") map$log_MSY <- factor(NA) else map$log_FMSY <- factor(NA)
    }
	obj2 <- MakeADFun(data = Assessment@info$data, parameters = params, map = map, inner.control = Assessment@info$inner.control,
                      random = Assessment@obj$env$random, DLL = "MSEtool", silent = TRUE)
    opt2 <- optimize_TMB_model(obj2, Assessment@info$control)[[1]]

    if(!is.character(opt2)) nll <- opt2$objective else nll <- NA
    return(nll)
  }
  nll <- vapply(1:nrow(profile_grid), profile_fn, numeric(1), Assessment = Assessment, params = params, map = map) - Assessment@opt$objective
  profile_grid$nll <- nll

  if(figure) {
    FMSY.MLE <- Assessment@FMSY
    MSY.MLE <- Assessment@MSY
    if(joint_profile) {
      z.mat <- acast(profile_grid, list("FMSY", "MSY"), value.var = "nll")
      contour(x = FMSY, y = MSY, z = z.mat, xlab = expression(F[MSY]), ylab = "MSY", nlevels = 20)
      points(FMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
    } else {
      if(profile_par == "FMSY") xlab <- expression(F[MSY]) else xlab <- "MSY"
      plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
    }

    if(save_figure) {
      Model <- Assessment@Model
      prepare_to_save_figure()

      create_png(file.path(plot.dir, "profile_likelihood.png"))
      if(joint_profile) {
        z.mat <- acast(profile_grid, list("FMSY", "MSY"), value.var = "nll")
        contour(x = FMSY, y = MSY, z = z.mat, xlab = expression(F[MSY]), ylab = "MSY", nlevels = 20)
        points(FMSY.MLE, MSY.MLE, col = "red", cex = 1.5, pch = 16)
		msg <- "Joint profile likelihood of FMSY and MSY. Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate."
      } else {
        if(profile_par == "FMSY") xlab <- expression(F[MSY]) else xlab <- "MSY"
        plot(getElement(profile_grid, profile_par), profile_grid$nll, xlab = xlab, ylab = "Change in neg. log-likeilhood value", typ = "o", pch = 16)
		msg <- paste0("Profile likelihood of ", profile_par, ". Numbers indicate change in negative log-likelihood relative to the minimum. Red point indicates maximum likelihood estimate.")
      }
      dev.off()

      profile.file.caption <- c("profile_likelihood.png", msg)

      html_report(plot.dir, model = "Surplus Production (State-Space)", captions = matrix(profile.file.caption, nrow = 1),
                  name = Assessment@Name, report_type = "Profile_Likelihood")
      browseURL(file.path(plot.dir, "Profile_Likelihood.html"))
    }
  }
  return(profile_grid)
}




#' @importFrom gplots rich.colors
retrospective_SP_SS <- function(Assessment, nyr, figure = TRUE) retrospective_SP(Assessment, nyr, figure, TRUE)

