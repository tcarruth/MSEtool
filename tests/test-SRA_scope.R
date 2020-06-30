
#### For data-limited situations - see notes below

library(MSEtool)
library(testthat)

OM <- testOM
OM@nsim = 2

Hist = runMSE(OM, Hist = TRUE)

bool <- c(TRUE, FALSE)
data_master <- expand.grid(cond = c("catch", "effort"), Catch_with_effort = bool,
                           Index = bool, ML = bool, CAL = bool, CAA = bool)

data_master$CAA[data_master$ML] <- data_master$CAL[data_master$ML] <- NA

#### Function to run SRA_scope in parallel
lapply_fn <- function(i, data_master, OM, Hist) {
  args <- list(OM = OM)
  args$data <- list()

  if(data_master$cond == "catch") { # Condition on catch
    args$condition <- "catch"
    args$data$Chist <- Hist@Data@Cat[1, ]
  } else {
    if(data_master$Catch_with_effort[i]) { # Conditioned on effort, but may include some catches
      args$data$Chist <- Hist@Data@Cat[1, ]
      args$data$Chist[c(1:40)] <- NA
    }
    args$condition <- "effort"
    args$data$Ehist <- Hist@TSdata$Find[1, ]
  }

  if(data_master$Index[i]) { # Only recent ten years index
    args$data$Index <- Hist@Data@Ind[1, ]
    args$data$Index[c(1:40)] <- NA
    args$data$I_sd <- c(rep(NA, 40), rep(0.3, 10))
    args$s_selectivity <- "B"
  }

  if(data_master$ML[i]) { # Only recent mean length
    args$data$ML <- Hist@Data@ML[1, ]
    args$data$ML[c(1:40)] <- NA
  } else if(data_master$CAA[i]) {

    args$data$CAA <- Hist@Data@CAA[1, , ]
    args$data$CAA[1:35, ] <- NA
  } else if(data_master$CAL[i]) {

    args$data$CAL <- Hist@Data@CAL[1, , ]
    args$data$CAL[1:35, ] <- NA

    args$data$length_bin <- Hist@Data@CAL_mids
  }
  args$mean_fit <- TRUE
  SRA <- do.call(SRA_scope, args)

  return(SRA)
}


DLMtool::setup(6)
sfExportAll()






#### Run SRA
res <- sfClusterApplyLB(1:nrow(data_master), lapply_fn, data_master = data_master, Hist = Hist, OM = OM)


#### Test whether plot function works
lapply_fn_plot <- function(x, res) {
  try(R.utils::withTimeout(plot(res[[x]], compare = FALSE, open_file = FALSE, filename = as.character(x)),
                           timeout = 5),
      silent = TRUE)
}

#sfLibrary(R.utils)
#sfExport(list = c("lapply_fn_plot", "res"))
test_plot <- lapply(1:nrow(data_master), lapply_fn_plot, res = res)
saveRDS(test_plot, file = "tests/test_plot.rds")

error_fn <- function(x) {
  length(grep("elapsed time limit", x[1])) > 0 | length(grep("C:/", x[1])) > 0
}



data_master$conv <- lapply(res, function(x) sum(x@conv)/length(x@conv))
data_master$figure <- vapply(test_plot, error_fn, numeric(1))
View(data_master)

