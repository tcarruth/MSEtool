
#' Remove observation error, process error, implementation error or future gradients in (time varying) parameters
#'
#' This function allows like-by-like comparison of results among DLMtool runMSE() and MSEtool multiMSE()
#'
#' @param MOM Object of class MOM. A Multi-OM object
#' @param obs Logical. Whether observation error should be removed
#' @param imp Logical. Whether implementation error should be removed
#' @param proc Logical. Whether process error should be removed (e.g. interannual variability in natural mortality rate)
#' @param grad Logical. Whether gradients (consistent temporal changes in parameters) should be removed
#' @param silent Logical. Should changes be printed to the console?
#' @author T.Carruthers and A. Hordyk
#' @keywords internal
#' @export
tinyErr_MOM<-function (MOM, obs = TRUE, imp = TRUE, proc = TRUE, grad = TRUE,
          silent = FALSE)
{
  if (!inherits(MOM, "MOM"))
    stop("Object must be class `MOM`", call. = FALSE)

  MOMout <- MOM

  np<-length(Stocks)
  nf<-length(Fleets[[1]])

  if (obs) {
    if (!silent)
      message("Removing all Observation Error")
    for(p in 1:np){for(f in 1:nf){
      MOMout@Obs[[p]][[f]]<-DLMtool::Perfect_Info
    }}

  }
  if (imp) {
    if (!silent)
      message("Removing all Implementation Error")
    for(p in 1:np){for(f in 1:nf){
      MOMout@Imps[[p]][[f]]<-DLMtool::Perfect_Imp
    }}
  }
  if (proc) {
    if (!silent)
      message("Removing all Process Error")

    vars <- c("cv", "sd", "Perr")

    # Stock P error
    nms <- c(slotNames("Stock"))
    ind <- unique(grep(paste(vars, collapse = "|"), nms,
                       value = FALSE))
    for(p in 1:np){

      for (X in seq_along(ind)) {
        n <- length(slot(MOMout@Stocks[[p]], nms[ind[X]]))
        if (n == 0)
          n <- 2
        slot(MOMout@Stocks[[p]], nms[ind[X]]) <- rep(0, n)
      }

    }

    # Fleet P error

    nms <- c(slotNames("Fleet"))
    ind <- unique(grep(paste(vars, collapse = "|"), nms,
                       value = FALSE))
    for(p in 1:np){
      for(f in 1:nf){
        for (X in seq_along(ind)) {
          n <- length(slot(MOMout@Fleets[[p]][[f]], nms[ind[X]]))
          if (n == 0)
            n <- 2
          slot(MOMout@Fleets[[p]][[f]], nms[ind[X]]) <- rep(0, n)
        }

      }
    }

  }

  if (grad) {
    if (!silent)
      message("Removing all Gradients")
      vars <- c("grad", "inc")

      # Stock grad
      nms <- c(slotNames("Stock"))
      ind <- unique(grep(paste(vars, collapse = "|"), nms,
                         value = FALSE))
      for(p in 1:np){

        for (X in seq_along(ind)) {
          n <- length(slot(MOMout@Stocks[[p]], nms[ind[X]]))
          if (n == 0)
            n <- 2
          slot(MOMout@Stocks[[p]], nms[ind[X]]) <- rep(0, n)
        }

      }

      nms <- c(slotNames("Fleet"))
      ind <- unique(grep(paste(vars, collapse = "|"), nms,
                         value = FALSE))
      for(p in 1:np){
        for(f in 1:nf){
          for (X in seq_along(ind)) {
            n <- length(slot(MOMout@Fleets[[p]][[f]], nms[ind[X]]))
            if (n == 0)
              n <- 2
            slot(MOMout@Fleets[[p]][[f]], nms[ind[X]]) <- rep(0, n)
          }

        }
      }


  }
  MOMout
}
