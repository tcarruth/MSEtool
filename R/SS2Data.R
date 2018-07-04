
#' Reads data Stock Synthesis file structure into an data object using package r4ss
#'
#' @description A function that uses the file location of a fitted SS3 model including input files to population
#' the various slots of an Data object.
#' @param SSdir A folder with Stock Synthesis input and output files in it
#' @param Source Reference to assessment documentation e.g. a url
#' @param Name The name for the Data object
#' @param Author Who did the assessment
#' @param comp_fleet A vector of indices over which to aggregate the composition (catch-at-length and catch-at-age) data.
#' \code{"all"} will aggregate across all fleets.
#' @param comp_partition Integer for selecting length/age observations that are retained, discarded, or both.
#' @param comp_gender Integer for selecting length/age observations that are female, male, or both.
#' @param index_fleet Integer for selecting the fleet of the index to put in the Data object.
#' @param depletion_basis Whether empirical or assessment based.
#' @param ... Arguments to pass to \link[r4ss]{SS_output}
#' @return An object of class Data.
#' @author T. Carruthers
#' @export
#' @seealso \link{SS2OM}
#' @importFrom r4ss SS_output
SS2Data <- function(SSdir, Source = "No source provided", Name = "", Author = "No author provided",
                    comp_fleet = "all", comp_partition = "all", comp_gender = "all", index_fleet = 1,
                    depletion_basis = c("empirical", "assess"), ...) {

  dots <- list(dir = SSdir, ...)
  if(!any(names(dots) == "covar")) dots$covar <- FALSE
  if(!any(names(dots) == "forecast")) dots$forecast <- FALSE
  if(!any(names(dots) == "ncols")) dots$ncols <- 1e3
  if(!any(names(dots) == "printstats")) dots$printstats <- FALSE
  if(!any(names(dots) == "verbose")) dots$verbose <- FALSE
  if(!any(names(dots) == "warn")) dots$warn <- FALSE

  message(paste("-- Using function SS_output of package r4ss version", packageVersion("r4ss"), "to extract data from SS file structure --"))
  replist <- do.call(SS_output, dots)
  message("-- End of r4ss operations --")

  season_as_years <- FALSE
  if(replist$nseasons == 1 && replist$seasduration < 1) {
    message(paste("Season-as-years detected in SS model. There is one season in the year with duration of", replist$seasduration, "year."))
    season_as_years <- TRUE
    nseas <- 1/replist$seasduration
    message("DLMtool operates on annual basis. Since the SS model is seasonal, we need to aggregate over seasons.")
  } else {
    nseas <- replist$nseasons
    if(nseas > 1) {
      message("DLMtool operating model is an annual model. Since the SS model is seasonal, we need to aggregate over seasons.")
    }
  }

  # Create Data object
  Data <- new("Data", stock = "MSE")
  if(is.null(Name)) {
    Data@Name <- SSdir
  } else Data@Name <- Name

  mainyrs <- replist$startyr:replist$endyr
  if(season_as_years) {
    nyears <- ceiling(length(mainyrs)/nseas)
    Data@Year <- 1:nyears
  } else {
    Data@Year <- mainyrs
  }

  ##### Life history
  #### Growth --------------------------------------
  growdat <- getGpars(replist)      # Age-specific parameters in endyr

  # Max age
  Data@MaxAge <- maxage <- ceiling(nrow(growdat)/ifelse(season_as_years, nseas, 1))

  seas1_aind_full <- expand.grid(nseas = 1:nseas, age = 1:maxage)
  seas1_aind <- which(seas1_aind_full$nseas == 1)

  GP <- replist$Growth_Parameters   # Some growth parameters (presumably in endyr)
  if(nrow(GP)>1){
    message(paste(nrow(GP),"different rows of growth parameters were reported by r4ss:"))
    print(GP)
    message("Only the first row of values will be used.\n")
  }

  #### Length at age --------------------------------------
  Data@vbLinf <- GP$Linf[1]
  Data@CV_vbLinf <- GP$CVmax[1]
  t0 <- GP$A_a_L0[1]
  #t0[t0 > 1] <- 0
  Data@vbt0 <- t0
  muK <- GP$K[1]
  if(muK <= 0) { #Estimate K from Len_age if K < 0 (e.g., age-varying K with negative deviations in K).
    get_K <- function(K, Lens, Linf, t0, ages) sum((Lens - (Linf * (1 - exp(-K * (ages - t0)))))^2)
    muK <- optimize(get_K, c(0, 2), Lens = Len_age, Linf = muLinf, t0 = t0, ages = 1:maxage)$minimum
  }
  Data@vbK <- muK

  message(paste0("Von Bertalanffy parameters: Linf = ", Data@vbLinf, ", K = ", muK, ", t0 = ", t0))

  #### Weight
  Data@wla <- GP$WtLen1[1]
  Data@wlb <- GP$WtLen2[1]

  message(paste0("Length-weight parameters: a = ", Data@wla, ", b = ", Data@wlb))

  #### Maturity --------------------------------------
  if(min(growdat$Len_Mat < 1)) {                    # Condition to check for length-based maturity
    Mat <- growdat$Len_Mat/max(growdat$Len_Mat)
  } else {                                          # Use age-based maturity
    Mat <- growdat$Age_Mat/max(growdat$Age_Mat)
  }
  if(season_as_years) Mat <- Mat[seas1_aind]

  # Currently using linear interpolation of mat vs len, is robust and very close to true logistic model predictions
  Len_age <- growdat$Len_Mid
  Data@L50 <- LinInterp(Mat, Len_age, 0.5+1e-6)
  Data@L95 <- LinInterp(Mat, Len_age, 0.95)

  message(paste0("Lengths at 50% and 95% maturity: ", paste(Data@L50, Data@L95, collapse = " ")))


  #### M --------------------------------------
  M <- growdat$M
  if(season_as_years) M <- M[seas1_aind]

  Data@Mort <- mean(M)
  if(length(unique(M)) > 1) {
    message(paste0("Age-varying M detected:\n", paste(M, collapse = " "), "."))
    message(paste0("Currently, DLMtool only supports a single value of M. Using mean to set Data@Mort = ", Data@Mort, "."))
  } else {
    message(paste0("Natural mortality M = ", Data@Mort))
  }

  #### Composition data -------------------------
  # (function argument to select comp or aggregate across all fleets/surveys)
  get_composition <- function(dbase) {
    dbase_ind <- match(dbase$Yr, mainyrs)
    dbase <- dbase[!is.na(dbase_ind), ]
    dbase$Obs2 <- dbase$Obs * dbase$N

    comp_list <- split(dbase, dbase$Fleet)
    comp_mat <- lapply(comp_list, acast, formula = list("Yr", "Bin"), fun.aggregate = sum, value.var = "Obs2", fill = 0)

    comp_ind <- match(comp_fleet, as.numeric(names(comp_mat)))
    comp_mat <- comp_mat[comp_ind]
    expand_matrix <- function(x) {
      res <- matrix(NA, nrow = length(mainyrs), ncol = maxage)
      res_ind <- match(as.numeric(rownames(x)), mainyrs)
      res[res_ind, ] <- x
      return(res)
    }
    comp_mat2 <- lapply(comp_mat, expand_matrix)
    comp_all <- do.call(rbind, comp_mat2)

    comp_res <- aggregate(comp_all, list(Yr = rep(1:length(mainyrs), length(comp_fleet))), sum, na.rm = TRUE)
    comp_res[comp_res == 0] <- NA

    return(comp_res[, -1])
  }

  # CAA
  if(nrow(agedbase) > 0) {
    if(is.character(comp_fleet) && comp_fleet == "all") comp_fleet <- unique(replist$agedbase$Fleet)
    if(any(!is.na(match(replist$agedbase$Fleet, comp_fleet)))) {
      CAA <- get_composition(replist$agedbase)
      Data@CAA <- array(as.matrix(CAA), c(1, nyears, maxage))
    } else {
      message(paste0("Could not find any age composition from Fleets: ", paste0(comp_fleet, collapse = " ")))
    }
    message(paste0("Collected age composition from Fleets: ", paste0(comp_fleet, collapse = " ")))
  }

  # CAL
  #if(nrow(lendbase) > 0) {
  #  if(is.character(comp_fleet) && comp_fleet == "all") comp_fleet <- unique(replist$lendbase$Fleet)
  #  if(any(!is.na(match(replist$lendbase$Fleet, comp_fleet)))) {
  #    CAL <- get_composition(replist$lendbase)
  #    Data@CAL <- array(as.matrix(CAL), c(1, nyears, ncol(CAL)))
  #    Data@CAL_bins <- replist$lbins # add one more length bin
  #  } else {
  #    message(paste0("Could not find any age composition from Fleets: ", paste0(comp_fleet, collapse = " ")))
  #  }
  #  message(paste0("Collected age composition from Fleets: ", paste0(comp_fleet, collapse = " ")))
  #}

  #### Catch -------------------------
  message("Adding total catch in weight across all fleets...")
  cat_yr_ind <- !is.na(match(replist$timeseries$Yr, mainyrs))
  ts <- replist$timeseries[cat_yr_ind, ]

  cat_col <- grepl("obs_cat", colnames(ts))
  cat <- ts[, cat_col]

  is_weight <- replist$catch_units[replist$IsFishFleet] == 1

  cat_weight <- cat[, is_weight]
  cat_numbers <- cat[, !is_weight]
  if(ncol(cat_numbers) > 0) {
    fleet_in_numbers <- which(replist$catch_units == 2)
    message(paste0("Catch in numbers was detected for Fleets: ", paste0(fleet_in_numbers, collapse = " ")))
    message("Will use estimated assessment output for catch in weight (dead biomass) for these Fleets.")

    cat_col2 <- grepl("dead\\(B", colnames(ts))
    cat2 <- ts[, cat_col2]
    cat_numbers <- cat2[, !is_weight]
  }
  total_catch <- aggregate(rowSums(cbind(cat_weight, cat_numbers)), list(Yr = ts$Yr), sum, na.rm = TRUE)

  tc_ind <- match(total_catch$Yr, mainyrs)
  total_catch_vec <- structure(total_catch$x[tc_ind], names = mainyrs)

  Data@Cat <- matrix(total_catch_vec, nrow = 1)
  message(paste0(sum(!is.na(tc_ind)), " years of catch in Data object."))
  message(paste0("Default CV of Catch = ", Data@CV_Cat))

  Data@AvC <- mean(total_catch)
  message(paste0("Mean catch Data@AvC = ", Data@AvC))

  #### Index -------------------------
  # (function argument to select index)
  cpue <- split(replist$cpue, replist$cpue$Fleet)
  cpue_ind <- match(index_fleet, as.numeric(names(cpue)))

  cpue_subset <- cpue[[cpue_ind]]

  cpue_ind2 <- match(cpue_subset$Yr, mainyrs)
  Ind <- structure(cpue_subset$Obs[cpue_ind2], names = mainyrs)
  if(is.na(Ind[length(Ind)])) warning("No index value in most recent year.")
  message(paste0(sum(!is.na(Ind)), " years of index values from Fleet ", index_fleet, " (", unique(cpue_subset$Fleet), ") in Data object."))
  Data@Ind <- matrix(Ind, nrow = 1)

  Data@CV_Ind <- mean(cpue_subset$SE/cpue_subset$Obs)
  message(paste0("CV of Index (as mean of the ratios of SE and observed values) = ", Data@CV_Ind))

  Data@Source <- paste0(Source,". Author: ", Author, ".")



  yind<-(1:nyears)*nseas
  yind2<-rep(1:nyears,each=nseas)



  cat<-rep(0,length(replist$timeseries$"obs_cat:_1"))
  for(i in 1:length(replist$timeseries)){
    if(grepl("obs_cat",names(replist$timeseries)[i]))cat<-cat+replist$timeseries[[i]]

  }
  if(nseas>1){

    ind<-rep(1:(length(cat)/nseas),nseas)
    cat2<-aggregate(cat,by=list(ind),sum)
    cat3<-cat2[1:length(yind2),2]
    cat4<-aggregate(cat3,by=list(yind2),sum)
  }


  Ind<-SSB
  Ind<-Ind/mean(Ind)
  Data@Ind <- matrix(Ind,nrow=1)
  rec<-replist$recruit$pred_recr
  rec<-rec[1:(nyears*nseas)]
  Recind<- aggregate(rec,by=list(yind2),mean)[,2]
  Recind<-Recind/mean(Recind)

  rec_ind <- match(mainyrs, replist$recruit$year)
  rec <- replist$recruit$pred_recr
  Data@Rec <- matrix(rec/mean(rec), nrow = 1)

  message("Relative recruitment strength (Data@Rec) obtained from assessment output.")

  Data@t <- nyears

  if(depletion == "assessment") {
    Data@Dep <- replist$current_depletion
    message(paste0("Based on assessment, depletion Data@Dep = ", Data@Dep))
  }
  Data@Dt<-Data@Dep<-Ind[nyears]/Ind[1]


  #### Reference points ----------------------
  Data@Cref <- replist$derived_quants$Value[replist$derived_quants$LABEL == "TotYield_MSY"]

  FMSY <- replist$derived_quants$Value[replist$derived_quants$LABEL == "Fstd_MSY"]
  Data@FMSY_M <- FMSY/Data@Mort

  Data@Bref <- replist$derived_quants$Value[replist$derived_quants$LABEL == "SSB_MSY"]
  SSB0 <- replist$derived_quants$Value[replist$derived_quants$LABEL == "SSB_Unfished"]
  Data@BMSY_B0 <- Data@Bref/SSB0

  Data@Iref<-Data@Ind[1]*Data@Bref/SSB[1]


  #### Steepness --------------------------------------
  steep <- replist$parameters[grepl("steep", rownames(replist$parameters)), ]
  if(nrow(steep) == 1) {
    Data@steep <- steep$Value
  } else {
    SR_ind <- match(mainyrs, replist$recruit$year)
    SSB <- replist$recruit$spawn_bio[SR_ind]
    rec <- replist$recruit$pred_recr[SR_ind]

    res <- try(SSB0 <- replist$derived_quants[replist$derived_quants$LABEL == "SPB_Virgin", 2], silent = TRUE)
    if(inherits(res, "try-error")) SSB0 <- SSB[1]

    res2 <- try(R0 <- replist$derived_quants[replist$derived_quants$LABEL == "Recr_Virgin", 2], silent = TRUE)
    if(inherits(res, "try-error")) {
      surv <- c(1, exp(-cumsum(M[1:(maxage-1)])))
      SpR0 <- sum(Wt_age * Mat * surv)
      R0 <- SSB0/SpR0
    } else {
      SpR0 <- SSB0/(R0 * ifelse(season_as_years, nseas, 1))
    }
    if(replist$SRRtype == 3 || replist$SRRtype == 6) SR <- "BH"
    if(replist$SRRtype == 2) SR <- "Ricker"
    Data@steep <- mean(SRopt(100, SSB, rec, SpR0, plot = FALSE, type = SR), na.rm = TRUE)
  }
  message(paste0("Stepness = ", Data@steep))

  SpAbun_ind <- match(replist$endyr+1, replist$recruit$year)
  Data@SpAbun <- replist$recruit$spawn_bio[SpAbun_ind]

  # Vulnerable biomass
  ts <- replist$timeseries[replist$timeseries$Yr == replist$endyr + 1, ]
  vb_ind <- grepl("sel\\(B", colnames(ts))

  vb <- ts[, vb_ind]
  Data@Abun <- sum(vb)


  ages<-growdat$Age
  cols<-match(ages,names(replist$Z_at_age))
  F_at_age=t(replist$Z_at_age[,cols]-replist$M_at_age[,cols])
  F_at_age[nrow(F_at_age),]<-F_at_age[nrow(F_at_age)-1,]# ad-hoc mirroring to deal with SS missing predicitons of F in terminal age
  Ftab<-cbind(expand.grid(1:dim(F_at_age)[1],1:dim(F_at_age)[2]),as.vector(F_at_age))

  if(nseas>1){
    sumF<-aggregate(Ftab[,3],by=list(aind[Ftab[,1]],Ftab[,2]),mean,na.rm=T)
    sumF<-aggregate(sumF[,3],by=list(sumF[,1],yind2[sumF[,2]]),sum,na.rm=T)
  }else{
    sumF<-Ftab
  }

  sumF<-sumF[sumF[,2]<nyears,] # generic solution: delete partial observation of final F predictions in seasonal model (last season of last year is NA)
  V <- array(0, dim = c(maxage, nyears))
  V[,1:(nyears-1)]<-sumF[,3] # for some reason SS doesn't predict F in final year


  Find<-apply(V,2,max,na.rm=T) # get apical F

  ind<-as.matrix(expand.grid(1:maxage,1:nyears))
  V[ind]<-V[ind]/Find[ind[,2]]

  # guess at length parameters # this is over ridden anyway
  ind<-((nyears-3)*nseas)+(0:((nseas*3)-1))
  muFage<-as.vector(apply(F_at_age[,ind],1,mean))
  Vuln<-muFage/max(muFage,na.rm=T)

  Data@LFC<-LinInterp(Vuln,Len_age,0.05,ascending=T,zeroint=T)                            # not used if V is in cpars
  Data@LFS<-Len_age[which.min((exp(Vuln)-exp(1.05))^2 * 1:length(Vuln))]  # not used if V is in cpars

  ages<-growdat$Age
  cols<-match(ages,names(replist$Z_at_age))


  cols<-match(ages,names(replist$catage))
  CAA=as.matrix(replist$catage[,cols])
  yr<-rep(replist$catage$Yr,dim(CAA)[2])
  age<-rep(ages,each=dim(CAA)[1])
  CAAagg<-aggregate(as.vector(CAA),by=list(yr,age),sum)
  CAAagg<-CAAagg[CAAagg[,2]!=0,]
  CAAagg<-CAAagg[CAAagg[,1]<=(nyears*nseas),]
  CAAagg<-CAAagg[CAAagg[,2]<=length(aind),]

  yind3<-yind2[CAAagg[,1]]
  aind3<-aind[CAAagg[,2]]

  CAA2<-aggregate(CAAagg[,3],by=list(yind3,aind3),sum)
  Data@CAA<-array(0,c(1,nyears,maxage))
  Data@CAA[as.matrix(cbind(rep(1,nrow(CAA2)),CAA2[,1:2]))]<-CAA2[,3]



  # CAL data
  Data@CAL_bins<-c(0,replist$lbins)

  Binno<-match(replist$lendbase$Bin,Data@CAL_bins)
  Yrno<-yind2[replist$lendbase$Yr]
  CALagg<-aggregate(replist$lendbase$Obs*replist$lendbase$N,by=list(Yrno,Binno),sum)
  Data@CAL<-array(0,c(1,nyears,length(Data@CAL_bins)))
  Data@CAL[as.matrix(cbind(rep(1,nrow(CALagg)),CALagg[,1:2]))]<-CALagg[,3]
  Data@CAL<-ceiling(Data@CAL)

  Data@ML<-matrix(apply(Data@CAL[1,,]*rep(Data@CAL_bins,each=nyears),1,mean),nrow=1)
  Data@ML[Data@ML==0]<-NA
  lcpos<-apply(Data@CAL[1,,],1,which.max)
  Data@Lc<-matrix(Data@CAL_bins[lcpos],nrow=1)
  Data@Lc[Data@Lc==0]<-NA

  Data@Lbar<-matrix(rep(NA,nyears),nrow=1)

  for(i in 1:nyears){
    if(!is.na(Data@Lc[i])){
      ind<-Data@CAL_bins>Data@Lc[i]
      Data@Lbar[1,i]<-mean(Data@CAL[1,i,ind]*Data@CAL_bins[ind])
    }
  }

  Data@Abun<-Data@SpAbun<-SSB[nyears]

  Data@vbLinf=GP$Linf[1]
  Data@CV_vbLinf=GP$CVmax[1]
  Data@vbK=GP$K[1]*nseas
  Data@vbt0=GP$A_a_L0
  Data@wla=GP$WtLen1[1]
  Data@wlb=GP$WtLen2[1]
  SSBpR<-SSB[1]/rec[1]
  hs<-mean(SRopt(100,SSB,rec,SSBpR,plot=F,type="BH"))
  Data@steep<-hs

  Data@Units<-""
  Data@Ref<-mean(replist$derived_quants[grepl("OFLCatch",replist$derived_quants$LABEL),2])
  Data@Ref_type = "OFL"
  Data@LHYear<-nyears
  Data@MPeff<-1
  Data@Log<-list(warning="This data was created by an alpha version of SS2Data and should not be trusted!")
  Data@Misc<-list(warning="This data was created by an alpha version of SS2Data and should not be trusted!")
  Data@Name<-"This data was created by an alpha version of SS2Data and should not be trusted!"

  Data@MPrec<-Data@Ref<-Data@Cat[1,nyears]
  #Data@Ref<-mean(replist$derived_quants[grepl("OFLCatch",replist$derived_quants$LABEL),2])
  Data@Ref_type = "Most recent annual catch"
  Data

}

