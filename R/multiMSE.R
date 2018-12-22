# =================================================================================================
# ====== multiMSE  - an R function for doing MSE evaluation for a multi-stock multi-fleet OM ======
# =================================================================================================

#' Run a multi-fleet multi-stock Management Strategy Evaluation
#'
#' A function that runs a Management Strategy Evaluation (closed-loop
#' simulation) for a specified operating model
#'
#'
#' @param MOM A multi-fleet multi-stock operating model (class 'MOM')
#' @param MPs A matrix of methods (nstock x nfleet) (character string) of class MP
#' @param CheckMPs Logical to indicate if \link{Can} function should be used to check
#' if MPs can be run.
#' @param timelimit Maximum time taken for a method to carry out 10 reps
#' (methods are ignored that take longer)
#' @param Hist Should model stop after historical simulations? Returns a list
#' containing all historical data
#' @param ntrials Maximum of times depletion and recruitment deviations are
#' resampled to optimize for depletion. After this the model stops if more than
#' percent of simulations are not close to the required depletion
#' @param fracD Maximum allowed proportion of simulations where depletion is not
#' close to sampled depletion from OM before model stops with error
#' @param CalcBlow Should low biomass be calculated where this is the spawning
#' biomass at which it takes HZN mean generation times of zero fishing to reach
#' Bfrac fraction of SSBMSY
#' @param HZN The number of mean generation times required to reach Bfrac SSBMSY
#' in the Blow calculation
#' @param Bfrac The target fraction of SSBMSY for calculating Blow
#' @param AnnualMSY Logical. Should MSY statistics be calculated for each projection year?
#' May differ from MSY statistics from last historical year if there are changes in productivity
#' @param silent Should messages be printed out to the console?
#' @param PPD Logical. Should posterior predicted data be included in the MSE object Misc slot?
#' @param parallel Logical. Should the MSE be run using parallel processing?
#' @param save_name Character. Optional name to save parallel MSE list
#' @param checks Logical. Run tests?
#' @param control control options for testing and debugging
#' @param maxsims Maximum number of simulations per packet
#' @param name Character string for name of saved MSE packets (if \code{savePack=TRUE})
#' and final MSE object. If none provided, it uses the first five letters from the \code{OM} name
#' @param unique Logical. Should the name be unique? Current date and time appended to name.
#' @param maxCrash Maximum number of consecutive crashes before the MSE stops
#' @param saveMSE Logical to indicate if final MSE object should be saved to current
#' working directory (this is probably a good idea)
#' @param savePack Logical to indicate if packets should be save to current working directory
#' @param ... Arguments to runMSE function
#' @return A hierarchical list (by stock then fleet) of objects of class MSE \linkS4class{MSE}
#' @author T. Carruthers and A. Hordyk
#' @export
multiMSE <- function(MOM, MPs = list(c("AvC","DCAC"),c("FMSYref","curE")),
                   CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=TRUE,
                   HZN=2, Bfrac=0.5, AnnualMSY=TRUE, silent=FALSE, PPD=FALSE, parallel=FALSE,
                   save_name=NULL, checks=FALSE, control=NULL) {

  # Known issues:
  # 1) maxage is the same across species so computationally inefficient
  # 2) initdist is no yet working as a cpars value (indeed I'm not sure where this comes into the OM specification at all!)

  if (class(MOM)!='MOM') stop("OM object is not of class '<OM'", call. = FALSE)
  if (Hist & parallel) {
    message("Sorry! Historical simulations currently can't use parallel.")
    parallel <- FALSE
  }

  # Set DLMenv to be empty. Currently updated by Assess models in MSEtool
  rm(list = ls(DLMenv), envir = DLMenv)

  # check if custom MP names already exist in DLMtool
  tt <-  suppressWarnings(try(lsf.str(envir=globalenv()), silent=TRUE))
  if (class(tt)!="try-error") {
    gl.funs <- as.vector(tt)
    pkg.funs <- as.vector(ls.str('package:DLMtool'))
    if ('package:MSEtool' %in% search()) pkg.funs <- c(pkg.funs, as.vector(ls.str('package:MSEtool')))
    if (length(gl.funs)>0) {
      gl.clss <- unlist(lapply(lapply(gl.funs, get), class))
      gl.MP <- gl.funs[gl.clss %in% 'MP']
      if (length(gl.MP)>0) {
        inc.gl <- gl.MP[gl.MP %in% MPs]
        if (length(inc.gl)>0) {
          dup.MPs <- inc.gl[inc.gl %in% pkg.funs]
          if (length(dup.MPs)>0) {
            stop("Custom MP names already in DLMtool: ", paste0(dup.MPs, " "), "\nRename Custom MPs")
          }
        }
      }
    }
  }

  # Check MPs
  MPvec<-unlist(MPs)
  if (!all(is.na(MPvec))) {
    for (mm in MPvec) {
      chkMP <- try(get(mm), silent=TRUE)
      if (class(chkMP) != 'MP') stop(mm, " is not a valid MP", call.=FALSE)
    }
  }

  if (parallel) {

    if(!snowfall::sfIsRunning()) {
      # stop("Parallel processing hasn't been initialized. Use 'setup'", call. = FALSE)
      message("Parallel processing hasn't been initialized. Calling 'setup()' now")
      setup()
    }

    if (all(is.na(MPs))) MPs <- avail("MP")
    # Export Custom MPs #
    cMPs <- MPs[!MPs %in% pkg.funs]
    if (length(cMPs)>0) snowfall::sfExport(list=cMPs)

    ncpu <- snowfall::sfCpus()

    if (OM@nsim<48) stop("nsim must be >=48")
    nits <- ceiling(OM@nsim/48)

    itsim <- rep(48,nits)

    if (nits < ncpu) {
      nits <- ncpu
      itsim <- rep(ceiling(OM@nsim/ncpu), ncpu)
    }
    cnt <- 1
    while(sum(itsim) != OM@nsim | any(itsim<2)) {
      diff <-  OM@nsim - sum(itsim)
      if (diff >0) {
        itsim[cnt] <- itsim[cnt]+1
      }
      if(diff < 0) {
        itsim[cnt] <- itsim[cnt]-1
      }
      cnt <- cnt+1
      if (cnt > length(itsim)) cnt <- 1
    }


    if (!silent) message("Running MSE in parallel on ", ncpu, ' processors')

    temp <- snowfall::sfClusterApplyLB(1:nits, run_parallel, itsim=itsim, OM=OM, MPs=MPs,
                                       CheckMPs=CheckMPs, timelimit=timelimit, Hist=Hist, ntrials=ntrials,
                                       fracD=fracD, CalcBlow=CalcBlow,
                                       HZN=HZN, Bfrac=Bfrac, AnnualMSY=AnnualMSY, silent=TRUE, PPD=PPD,
                                       control=control, parallel=parallel)
    #assign_DLMenv() # grabs objects from DLMenv in cores, then merges and assigns to 'home' environment

    if (!is.null(save_name) && is.character(save_name)) saveRDS(temp, paste0(save_name, '.rdata'))

    MSE1 <- joinMSE(temp)
    if (class(MSE1) == "MSE") {
      if (!silent) message("MSE completed")
    } else {
      warning("MSE completed but could not join MSE objects. Re-run with `save_name ='MyName'` to debug")
    }
  }

  if (!parallel) {
    if (MOM@nsim > 48 & !silent & !Hist) message("Suggest using 'parallel = TRUE' for large number of simulations")
    MSE1 <- multiMSE_int(MOM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow,
                       HZN, Bfrac, AnnualMSY, silent, PPD, checks=checks, control=control)

  }

  return(MSE1)

}

multiMSE_int <- function(MOM, MPs=list(c("AvC","DCAC"),c("FMSYref","curE")),
                       CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=TRUE,
                       HZN=2, Bfrac=0.5, AnnualMSY=TRUE, silent=FALSE, PPD=FALSE, checks=FALSE,
                       control=NULL, parallel=FALSE) {

  # Priority list
  # Multi-fleet popdynCPP needed to properly account for fleet-specific MPAs, spatial targetting etc.

  # Known issues
  # 1) SSN[SPAYR] <- Nfrac[SA] * StockPars[[p]]$R0[S] * StockPars[[p]]$initdist[SAR]  this is Pinitdist[SR] in runMSE but I have no idea where Pinitdist comes from
  # 2) currently assume unfished vulnerability is equally weighted among fleets (first V calculation)
  # 3) Historical MPAs are calculated by weighting fleet specific MPA by their historical effort / q (ie apical F) - not exactly ideal as it should be by vulnerable biomass really
  # 4) Spatial targetting is currently calculated as the fleet-specific spatial targetting weighted by total F by fleet - not perfect at all!
  # 5) MICE mode: MSY calcs are fully dynamic and by year - we do not solve the long-term equilibrium MSY for MICE models
  # 6) No check for MOM correct object formatting yet
  # 7) MICE mode: SSB0 is constant (does not change with M for example) and based on long-term ecosystem average as specified in StockPars[[p]]$SSB0

  # Needs checking
  # 1) The vulnerability in the fleets does not max to 1
  # 2) When calculating aggregate retention in the multi_q_estimation function, how should max retention at age be calculated (what is the max value)?

  if (class(MOM) != "MOM") stop("You must specify a valid operating model of class MOM (multi operating model)")
  Misc<-new('list') # Blank miscellaneous slot created
  if("seed"%in%slotNames(MOM)) set.seed(MOM@seed) # set seed for reproducibility

  #OM <- updateMSE(OM)
  tiny <- 1e-15  # define tiny variable
  nsim<-MOM@nsim
  nyears <- MOM@Fleets[[1]][[1]]@nyears  # number of historical years
  proyears<-MOM@proyears
  allyears<-nyears+proyears
  Stocks<-MOM@Stocks
  Fleets<-MOM@Fleets
  Rel<-MOM@Rel
  CatchFrac<-MOM@CatchFrac
  np<-length(Stocks)
  nf<-length(Fleets[[1]])
  maxF<-MOM@maxF

  #MOM <- ChkObj(MOM) # Check that all required slots in OM object contain values
  if (proyears < 2) stop('OM@proyears must be > 1', call.=FALSE)
  ### Sampling OM parameters ###
  if(!silent) message("Loading operating model")

  # --- Sample custom parameters ----

  SampCpars<-list() # empty list
  # custom parameters exist - sample and write to list

  for(p in 1:np){
    SampCpars[[p]]<-list()
    for(f in 1:nf){
      if(length(cpars[[p]][[f]])>0){
        message(paste(Stocks[[p]]@Name," - ",Fleets[[p]][[f]]@Name))
        ncparsim<-cparscheck(cpars[[p]][[f]])   # check each list object has the same length and if not stop and error report
        SampCpars[[p]][[f]] <- SampleCpars(cpars[[p]][[f]], nsim, msg=!silent)
      }
    }
  }

  # All stocks and sampled parameters must have compatible array sizes (maxage)
  maxage_s<-unique(SIL(MOM@Stocks,"maxage"))
  if(length(maxage_s)>1)message(paste("Stocks of varying maximum ages have been specified, all simulations will run to",max(maxage_s),"ages"))
  maxage<-max(maxage_s)
  for(p in 1:np)MOM@Stocks[[p]]@maxage<-maxage

  # --- Sample Stock Parameters ----
  StockPars<-FleetPars<-ObsPars<-ImpPars<-new('list')
  for(p in 1:np){
    StockPars[[p]] <- SampleStockPars(MOM@Stocks[[p]], nsim, nyears, proyears,  SampCpars[[p]][[1]], msg=!silent)
  }

  # --- Sample Fleet Parameters ----
  for(p in 1:np){
    FleetPars[[p]]<-ObsPars[[p]]<-ImpPars[[p]]<-list()
    for(f in 1:nf){
      FleetPars[[p]][[f]] <- SampleFleetPars(MOM@Fleets[[p]][[f]], Stock=StockPars[[p]], nsim, nyears, proyears, cpars=SampCpars[[p]][[f]])
      ObsPars[[p]][[f]] <- SampleObsPars(MOM@Obs[[p]][[f]], nsim, cpars=SampCpars[[p]][[f]])
      ImpPars[[p]][[f]] <- SampleImpPars(MOM@Imps[[p]][[f]],nsim, cpars=SampCpars[[p]][[f]])

      FleetPars[[p]][[f]]$Find<-FleetPars[[p]][[f]]$Find/(apply(FleetPars[[p]][[f]]$Find,1,mean))
    }
  }

  # Assign Stock pars to function environment
  #for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])

  # Assign Fleet pars to function environment
  #for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])

  # Assign Obs pars to function environment
  # for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])

  # Assign Imp pars to function environment
  # for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])

  # Check for compatibility of dimensions among stocks

  nareas_s<-NIL(StockPars,"nareas")
  if(length(unique(nareas_s))!=1)stop("Stocks must have the same specified number of areas - check cpars$mov for each stock object")
   nareas<-nareas_s[1]

  ### End of sampling OM parameters ###
  # nsim, ns, nf, maxage, nyears, nareas
  N <- Biomass <- Z<- VBiomass<- SSN <- SSB <- array(NA, dim = c(nsim, np, maxage, nyears, nareas))
  VF<-FretA<-array(NA, dim = c(nsim, np, nf, maxage, allyears))
  VBF<-FM <- FMret <- array(NA, dim = c(nsim, np, nf, maxage, nyears, nareas))  # stock numbers array
  SPR <- array(NA, dim = c(nsim, np, maxage, nyears)) # store the Spawning Potential Ratio
  MPA<-array(1,c(np,nf, nyears+proyears,nareas))
  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array

  # Survival array with M-at-age

  for(p in 1:np){

    surv <- matrix(1, nsim, maxage)
    surv[, 2:maxage] <- t(exp(-apply(StockPars[[p]]$M_ageArray[,,1], 1, cumsum)))[, 1:(maxage-1)]  # Survival array

    Nfrac <- surv * StockPars[[p]]$Mat_age[,,1]  # predicted Numbers of mature ages in first year

    SPAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, p, 1:nsim)[5:1])  # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
    SAY <- SPAYR[, c(1,3,4)]
    SAR <- SPAYR[, c(1,3,5)]
    SA <- Sa <- SPAYR[, c(1,3)]
    SR <- SPAYR[, c(1,5)]
    S <- SPAYR[, 1]
    SY <- SPAYR[, c(1, 4)]
    Sa[,2]<-maxage-Sa[,2]+1 # This is the process error index for initial year

    # <<<< Previously an if statement for existence of initdist, see runMSE.R of DLMtool/R>>>>>>>>

    R0a <- matrix(StockPars[[p]]$R0, nrow=nsim, ncol=nareas, byrow=FALSE) * StockPars[[p]]$initdist[,1,]  # !!!! INITDIST OF AGE 1. Unfished recruitment by area

    SSN[SPAYR] <- Nfrac[SA] * StockPars[[p]]$R0[S] * StockPars[[p]]$initdist[SAR]  # Calculate initial spawning stock numbers
    N[SPAYR] <- StockPars[[p]]$R0[S] * surv[SA] * StockPars[[p]]$initdist[SAR]  # Calculate initial stock numbers
    Neq <- N
    Biomass[SPAYR] <- N[SPAYR] * StockPars[[p]]$Wt_age[SAY]  # Calculate initial stock biomass
    SSB[SPAYR] <- SSN[SPAYR] * StockPars[[p]]$Wt_age[SAY]    # Calculate spawning stock biomass

    Vraw<-array(NIL(listy=FleetPars[[p]],namey="V"),c(nsim,maxage,allyears,nf))
    Vind<-as.matrix(expand.grid(1:nsim,p,1:nf,1:maxage,1:allyears))
    VF[Vind]<-Vraw[Vind[,c(1,4,5,3)]]
    if(nf==1){
      V<-VF[,p,1,,] #<-SOL(FleetPars[[p]],"V")
    }else{
      V<-apply(VF[,p,,,],c(1,3,4),sum)
    }
    V<-nlz(V,c(1,3),"max") # currently assume unfished vulnerability is equally weighted among fleets
    VBiomass[SPAYR] <- Biomass[SPAYR] * V[SAY]  # Calculate vunerable biomass

    Fretraw<-array(NIL(listy=FleetPars[[p]],namey="retA"),c(nsim,maxage,allyears,nf))
    FretA[Vind]<-Fretraw[Vind[,c(1,4,5,3)]]

    if (nsim > 1) {
      SSN0 <- apply(SSN[,p , , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
      SSB0 <- apply(SSB[,p , , 1, ], 1, sum)  # Calculate unfished spawning stock biomass
      SSBpR <- matrix(SSB0/StockPars[[p]]$R0, nrow=nsim, ncol=nareas)  # Spawning stock biomass per recruit
      SSB0a <- apply(SSB[,p, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
      B0 <- apply(Biomass[,p, , 1, ], 1, sum)
      N0 <- apply(N[,p, , 1, ], 1, sum)
    } else {
      SSN0 <- apply(SSN[,p, , 1, ], 2, sum)  # Calculate unfished spawning stock numbers
      SSB0 <-  sum(SSB[,p, , 1, ])  # Calculate unfished spawning stock biomass
      SSBpR <- SSB0/R0  # Spawning stock biomass per recruit
      SSB0a <- apply(SSB[,p, , 1, ], 2, sum)  # Calculate unfished spawning stock numbers
      B0 <- apply(Biomass[,p, , 1, ], 2, sum)
      N0 <- apply(N[,p, , 1, ], 2, sum)
    }

    bR <- matrix(log(5 * StockPars[[p]]$hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
    aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params

    #  --- Non-equilibrium calcs ----
    SSN[SPAYR] <- Nfrac[SA] * StockPars[[p]]$R0[S] * StockPars[[p]]$initdist[SAR]*StockPars[[p]]$Perr_y[Sa]  # Calculate initial spawning stock numbers
    N[SPAYR] <- StockPars[[p]]$R0[S] * surv[SA] * StockPars[[p]]$initdist[SAR]*StockPars[[p]]$Perr_y[Sa]  # Calculate initial stock numbers

    Biomass[SPAYR] <- N[SPAYR] * StockPars[[p]]$Wt_age[SAY]  # Calculate initial stock biomass
    SSB[SPAYR] <- SSN[SPAYR] * StockPars[[p]]$Wt_age[SAY]    # Calculate spawning stock biomass
    VBiomass[SPAYR] <- Biomass[SPAYR] * V[SAY]  # Calculate vunerable biomass

    # Assign stock parameters to StockPars object

    StockPars[[p]]$SSBpR <-SSBpR
    StockPars[[p]]$aR <-aR
    StockPars[[p]]$bR<-bR
    StockPars[[p]]$SSB0<-SSB0
    StockPars[[p]]$R0a <-R0a
    StockPars[[p]]$surv<-surv


    for(f in 1:nf)FleetPars[[p]][[f]]$V<-VF[,p,f,,]

    # --- Historical Spatial closures ----

    for(f in 1:nf){

      if (all(!is.na(Fleets[[p]][[f]]@MPA)) && sum(Fleets[[p]][[f]]@MPA) != 0) {

        yrindex <- Fleets[[p]][[f]]@MPA[,1]
        if (max(yrindex)>nyears) stop("Invalid year index for spatial closures: must be <= nyears")
        if (min(yrindex)<1) stop("Invalid year index for spatial closures: must be > 1")
        if (ncol(Fleets[[p]][[f]]@MPA)-1 != nareas) stop("OM@MPA must be nareas + 1")

        for (xx in seq_along(yrindex)) {
          MPA[p,f,yrindex[xx]:nrow(MPA),] <- matrix(Fleets[[p]][[f]]@MPA[xx, 2:ncol(Fleets[[p]][[f]]@MPA)], nrow=length(yrindex[xx]:nrow(MPA)),ncol=nareas, byrow = TRUE)
        }

      }
    }

  } # end of np





  # --- Optimize catchability (q) to fit depletion ----
  #if ('unfished' %in% names(control) && control$unfished) {
  # if(!silent) message("Simulating unfished historical period")
  #  Hist <- TRUE
  # CalcBlow <- FALSE
  #  qs <- array(0,dim(nsim,nf)) # no fishing
  #} else {

  if(!silent) message("Optimizing for user-specified depletion")  # Print a progress update

  bounds <- c(0.0001, 15) # q bounds for optimizer
  out<-sfLapply(1:nsim,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF, MPA,CatchFrac, bounds= bounds,Rel)

  #out<-lapply(1:nsim,getq_multi_MICE, StockPars, FleetPars, nareas,maxage, pyears, N, maxF=MOM@maxF, MPA, CatchFrac, bounds)
  #return and reallocate
  qs<-t(matrix(NIL(out,"qtot"),nrow=np))
  qout<-array(NIL(out,"qfrac"),c(np,nf,nsim))
  qfrac<-array(NA,c(nsim,np,nf))
  qind<-TEG(dim(qfrac))
  qfrac[qind]<-qout[qind[,c(2,3,1)]]

    #out<-lapply(1:sim,getq_multi, D=StockPars[[p]]$D, SSB0=SSB0, nareas=nareas, maxage=maxage, Np=N[,p,,,], pyears=nyears, M_ageArray=StockPars[[p]]$M_ageArray, Mat_age=StockPars[[p]]$Mat_age,
    #                Asize=StockPars[[p]]$Asize, Wt_age=StockPars[[p]]$Wt_age,
    #                FleetP=FleetPars[[p]], # instead of V, retA, Find, Spat_targ,
    #                Perr=StockPars[[p]]$Perr_y, mov=StockPars[[p]]$mov, SRrel=StockPars[[p]]$SRrel, hs=StockPars[[p]]$hs, R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR,
    #                bounds = c(1e-05, 15), maxF=MOM@maxF, MPAc=MPA[p,,,], CFp=CatchFrac[[p]], useCPP=TRUE)

    #qs<-NIL(out,"qtot")
    #qfrac<-t(matrix(NIL(out,"qfrac"),nrow=nf))


  # --- Check that q optimizer has converged ----
  #qs[3,2]<-qs[2,1]<- -10000
  LimBound <- c(1.1, 0.9)*range(bounds)  # bounds for q (catchability). Flag if bounded optimizer hits the bounds
  probQ <- which(apply(qs > max(LimBound) | qs < min(LimBound),1,sum)>0)
  Nprob <- length(probQ)

  # If q has hit bound, re-sample depletion and try again. Tries 'ntrials' times and then alerts user
  if (length(probQ) > 0) {

    Err <- TRUE
    if(!silent) message(Nprob,' simulations have final biomass that is not close to sampled depletion')
    if(!silent) message('Re-sampling depletion, recruitment error, and fishing effort')

    count <- 0
    MOM2 <- MOM

    while (Err & count < ntrials) {
      # Re-sample Stock Parameters

      Nprob <- length(probQ)
      MOM2@nsim <- Nprob

      SampCpars2<-list()
      for(f in 1:nf){
        if(length(cpars[[p]][[f]])>0){
          #message(paste(Stocks[[p]]@Name," - ",Fleets[[p]][[f]]@Name))
          ncparsim<-cparscheck(cpars[[p]][[f]])   # check each list object has the same length and if not stop and error report
          SampCpars2[[f]] <- SampleCpars(cpars[[p]][[f]], Nprob, msg=!silent)
        }
      }

      #StockPars[[p]] <- SampleStockPars(MOM@Stocks[[p]], nsim, nyears, proyears,  SampCpars[[p]][[f]], msg=!silent)
      ResampStockPars <- SampleStockPars(MOM2@Stocks[[p]], nsim=Nprob,nyears=nyears,proyears=proyears,cpars=SampCpars2[[1]], msg=FALSE)
      #ResampStockPars$CAL_bins <- StockPars$CAL_bins
      #ResampStockPars$CAL_binsmid <- StockPars$CAL_binsmid

      # Re-sample depletion
      StockPars[[p]]$D[probQ] <- ResampStockPars$D

      # Re-sample recruitment deviations
      StockPars[[p]]$procsd[probQ] <- ResampStockPars$procsd
      StockPars[[p]]$AC[probQ] <- ResampStockPars$AC
      StockPars[[p]]$Perr_y[probQ,] <- ResampStockPars$Perr_y
      StockPars[[p]]$hs[probQ] <- ResampStockPars$hs

      # Re-sample historical fishing effort
      ResampFleetPars<-list()

      for(f in 1:nf){
        ResampFleetPars <- SampleFleetPars(MOM2@Fleets[[p]][[f]], Stock=ResampStockPars, nsim=Nprob, nyears=nyears, proyears=proyears, cpars=SampCpars2[[f]])
        FleetPars[[p]][[f]]$Esd[probQ] <- ResampFleetPars$Esd
        FleetPars[[p]][[f]]$Find[probQ, ] <- ResampFleetPars$Find
        FleetPars[[p]][[f]]$dFfinal[probQ] <- ResampFleetPars$dFfinal
      }


      out2<-sfLapply(probQ,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF, MPA,CatchFrac, bounds= bounds,Rel)

      qs2<-t(matrix(NIL(out2,"qtot"),nrow=np))
      qout2<-array(NIL(out2,"qfrac"),c(np,nf,nsim))
      qfrac2<-array(NA,c(Nprob,np,nf))
      qind2<-TEG(dim(qfrac2))
      qfrac2[qind2]<-qout2[qind2[,c(2,3,1)]]

      #qs2<-NIL(out2,"qtot")
      #qfrac2<-t(matrix(NIL(out2,"qfrac"),nrow=nf))
      qfrac[probQ,,]<-qfrac2
      qs[probQ,]<-qs2

      #if(nf>1){
        # Optimize for q

       # out2<-lapply(probQ,getq_multi, D=StockPars[[p]]$D, SSB0=SSB0, nareas=nareas, maxage=maxage, Np=N[,p,,,], pyears=nyears, M_ageArray=StockPars[[p]]$M_ageArray, Mat_age=StockPars[[p]]$Mat_age,
         #           Asize=StockPars[[p]]$Asize, Wt_age=StockPars[[p]]$Wt_age, FleetP=FleetPars[[p]], Perr=StockPars[[p]]$Perr_y, mov=StockPars[[p]]$mov, SRrel=StockPars[[p]]$SRrel, hs=StockPars[[p]]$hs, R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR,
        #            bounds = bounds, maxF=MOM@maxF, MPAc=MPA[p,,,], CFp=CatchFrac[[p]], useCPP=TRUE)

        #qs2<-NIL(out2,"qtot")
        #qfrac2<-t(matrix(NIL(out2,"qfrac"),nrow=nf))
        #qfrac[probQ,]<-qfrac2
        #qs[probQ]<-qs2

      #}else{

       # qs2 <- sapply(probQ, getq3, D=StockPars[[p]]$D, SSB0=SSB0, nareas=nareas, maxage=maxage, N=N[,p,,,], pyears=nyears,
        #             M_ageArray=StockPars[[p]]$M_ageArray, Mat_age=StockPars[[p]]$Mat_age, Asize=StockPars[[p]]$Asize,
         #            Wt_age=StockPars[[p]]$Wt_age, V=V, retA=FleetPars[[p]][[1]]$retA, Perr=StockPars[[p]]$Perr_y, mov=StockPars[[p]]$mov,
          #           SRrel=StockPars[[p]]$SRrel, Find=FleetPars[[p]][[1]]$Find,
           #          Spat_targ=FleetPars[[p]][[1]]$Spat_targ, hs=StockPars[[p]]$hs, R0a=R0a, SSBpR=SSBpR, aR=aR, bR=bR,
            #         bounds=bounds, MPA=MPA[p,1,,], maxF=MOM@maxF) # find the q that gives current stock depletion

        #qs[probQ]<-qs2
        # need a qfract thing in here

      #}

      probQ <- which(apply(qs > max(LimBound) | qs < min(LimBound),1,sum)>0)
      count <- count + 1
      if (length(probQ) == 0) Err <- FALSE

    }
    if (Err) { # still a problem

      tooLow <- length(which(qs > max(LimBound)))
      tooHigh <- length(which(qs < min(LimBound)))
      prErr <- length(probQ)/nsim
      if (prErr > fracD & length(probQ) >= 1) {
        if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
        if(!silent) message("More than ", fracD*100, "% of simulations can't get to the specified level of depletion with these Operating Model parameters")
        stop("Change OM@seed and try again for a complete new sample, modify the input parameters, or increase ntrials")
      } else {
        if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
        if(!silent) message("More than ", 100-fracD*100, "% simulations can get to the sampled depletion.\nContinuing")
      }
    }

    for(f in 1:nf) FleetPars[[p]][[ff]]$qs<-qs*qfrac[,f]

  }  # end of q estimation steps

  if(!silent) message("Calculating historical stock and fishing dynamics")  # Print a progress update

  # Got to here!


  # need retA, Effind, Vuln, Spat_targ


  histYrs <- sapply(1:nsim, HistMulti, FleetP=FleetPars[[p]], StockP=StockPars[[p]],
                    maxage=maxage,allyears=allyears, nyears=nyears, nf=nf,
                    MPAc=MPA[p,,,], qfrac=qfrac, qs=qs, maxF=maxF)


  N[,p,,,] <- aperm(array(as.numeric(unlist(histYrs[1,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
  Biomass[,p,,,] <- aperm(array(as.numeric(unlist(histYrs[2,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
  SSN[,p,,,] <- aperm(array(as.numeric(unlist(histYrs[3,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
  SSB[,p,,,] <- aperm(array(as.numeric(unlist(histYrs[4,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))
  VBiomass[,p,,,] <- aperm(array(as.numeric(unlist(histYrs[5,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))

  for(f in 1:nf){
    VBFind<-as.matrix(expand.grid(1:nsim,p,f,1:maxage,1:nyears,1:nareas))
    VBF[VBFind]<- Biomass[VBFind[,c(1:2,4:6)]] * FleetPars[[p]][[f]]$V[VBFind[,c(1,4,5)]]
  }
  #matplot(t(VBF[1,1,,,10,1]))
  #matplot(t(VBF[1:3,1,1,,10,1]))

  Effind<-array(NIL(FleetPars[[p]],"Find"),dim=c(nsim,nyears,nf))
  FMind<-as.matrix(expand.grid(1:nsim,p,1:nf,1:maxage,1:nyears,1:nareas))
  FM[FMind]<-Effind[FMind[,c(1,5,3)]]*qs[FMind[,1]]*qfrac[FMind[,c(1,3)]]*VF[FMind[,1:5]]
  FMret[FMind] <-  FM[FMind]*FretA[FMind[,1:5]]
  Z[,p,,,]<-aperm(array(as.numeric(unlist(histYrs[8,], use.names=FALSE)), dim=c(maxage, nyears, nareas, nsim)), c(4,1,2,3))

  if (nsim > 1) Depletion <- apply(SSB[,p,,nyears,],1,sum)/SSB0#^betas
  if (nsim == 1) Depletion <- sum(SSB[,p,,nyears,])/SSB0 #^betas


  # # apply hyperstability / hyperdepletion

  # Check that depletion is correct
  if (checks) {
    if (prod(round(StockPars[[p]]$D, 2)/ round(Depletion,2)) != 1) {
      print(cbind(round(D,4), round(Depletion,4)))
      warning("Possible problem in depletion calculations")
    }
  }

  # end of stocks
  #} # end of if not unfished

  if(!silent) message("Calculating MSY reference points")  # Print a progress update

  for(p in 1:np){
  # --- Calculate MSY references ----

    V<-apply(VF[,p,,,],c(1,3,4),sum) #<-SOL(FleetPars[[p]],"V")
    V<-nlz(V,c(1,3),"max")
    MSYrefs <- sapply(1:nsim, optMSY_eq, StockPars[[p]]$M_ageArray, StockPars[[p]]$Wt_age, StockPars[[p]]$Mat_age,
                      V=V, maxage=maxage,
                      R0=StockPars[[p]]$R0, SRrel=StockPars[[p]]$SRrel, hs=StockPars[[p]]$hs, yr=nyears)

    StockPars[[p]]$MSY <- MSYrefs[1, ]  # record the MSY results (Vulnerable)
    StockPars[[p]]$FMSY <- MSYrefs[2, ]  # instantaneous FMSY (Vulnerable)
    StockPars[[p]]$SSBMSY <- MSYrefs[3, ]  # Spawning Stock Biomass at MSY
    StockPars[[p]]$SSBMSY_SSB0 <- MSYrefs[4, ] # SSBMSY relative to unfished (SSB)
    StockPars[[p]]$BMSY_B0 <- MSYrefs[5, ] # Biomass relative to unfished (B0)
    StockPars[[p]]$BMSY <- MSYrefs[6,] # total biomass at MSY
    StockPars[[p]]$VBMSY <- (MSYrefs[1, ]/(1 - exp(-MSYrefs[2, ])))  # Biomass at MSY (Vulnerable)
    StockPars[[p]]$UMSY <- MSYrefs[1, ]/StockPars[[p]]$VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
    StockPars[[p]]$FMSY_M <- StockPars[[p]]$FMSY/StockPars[[p]]$M  # ratio of true FMSY to natural mortality rate M

    if (checks) {
      Btemp <- apply(SSB, c(1,3), sum)
      x <- Btemp[,nyears]/SSBMSY
      y <-D/SSBMSY_SSB0
      plot(x,y, xlim=c(0,max(x)), ylim=c(0,max(y)), xlab="SSB/SSBMSY", ylab="D/SSBMSY_SSB0")
      lines(c(-10,10),c(-10,10))
    }

  } # end of stocks

  # --- Code for deriving low biomass ----
  # (SSB where it takes MGThorizon x MGT to reach Bfrac of BMSY)

  Znow<-apply(Z[,,nyears,]*N[,,nyears,],1:2,sum)/apply(N[,,nyears,],1:2,sum)
  MGTsurv<-t(exp(-apply(Znow,1,cumsum)))
  MGT<-apply(Agearray*(Mat_age[,,nyears]*MGTsurv),1,sum)/apply(Mat_age[,,nyears]*MGTsurv,1,sum)

  if(CalcBlow){
    if(!silent) message("Calculating B-low reference points")              # Print a progress update

    MGThorizon<-floor(HZN*MGT)

    Blow <- sapply(1:nsim,getBlow, N, Asize, MSYrefs[3,],SSBpR, MPA, SSB0, nareas, retA, MGThorizon,
                   Find,Perr_y,M_ageArray,hs,Mat_age, Wt_age,R0a,V,nyears,maxage,mov,
                   Spat_targ,SRrel,aR,bR,Bfrac, maxF)
  }else{
    Blow<-rep(NA,nsim)
  }


  # --- Calculate Reference Yield ----
  if(!silent) message("Calculating reference yield - best fixed F strategy")  # Print a progress update

  RefY <- sapply(1:nsim, getFref3, Asize, nareas, maxage, N=N[,,nyears,, drop=FALSE], pyears=proyears,
                 M_ageArray=M_ageArray[,,(nyears):(nyears+proyears)], Mat_age[,,(nyears):(nyears+proyears)],
                 Wt_age=Wt_age[,,nyears:(nyears+proyears)],
                 V=retA[, , (nyears + 1):(nyears + proyears), drop=FALSE],
                 retA=retA[, , (nyears + 1):(nyears + proyears), drop=FALSE],
                 Perr=Perr_y[,(nyears):(nyears+maxage+proyears-1)], mov, SRrel, Find,
                 Spat_targ, hs, R0a, SSBpR, aR, bR, MPA=MPA, maxF=maxF, SSB0=SSB0)

  # --- Calculate catch-at-age ----
  CN <- apply(N * (1 - exp(-Z)) * (FM/Z), c(1, 3, 2), sum)  # Catch in numbers (removed from population)
  CN[is.na(CN)] <- 0
  CB <- Biomass * (1 - exp(-Z)) * (FM/Z)  # Catch in biomass (removed from population)

  # --- Calculate retained-at-age ----
  Cret <- apply(N * (1 - exp(-Z)) * (FMret/Z), c(1, 3, 2), sum)  # Retained catch in numbers
  Cret[is.na(Cret)] <- 0
  CBret <- Biomass * (1 - exp(-Z)) * (FMret/Z)  # Retained catch in biomass

  # --- Calculate dead discarded-at-age ----
  Cdisc <- CN - Cret # discarded numbers
  CBdisc <- CB - CBret # discarded biomass

  # --- Simulate observed catch ----

  if (length(control$Cbias_yr) ==0) {
    Cbiasa <- array(Cbias, c(nsim, nyears + proyears))  # Bias array
  } else {
    Cbiasa <- matrix(1, nsim, nyears+proyears)
    Cbiasa[,control$yrs] <- control$Cbias_yr
  }
  Cerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Csd, (nyears + proyears))),
                       sdconv(1, rep(Csd, nyears + proyears))), c(nsim, nyears + proyears))  # composite of bias and observation error
  # Cobs <- Cbiasa[, 1:nyears] * Cerr[, 1:nyears] * apply(CB, c(1, 3), sum)  # Simulated observed catch (biomass)
  Cobs <- Cbiasa[, 1:nyears] * Cerr[, 1:nyears] * apply(CBret, c(1, 3), sum)  # Simulated observed retained catch (biomass)

  # --- Simulate observed catch-at-age ----

  # generate CAA from retained catch-at-age
  CAA <- array(NA, dim = c(nsim, nyears, maxage))  # Catch  at age array
  # cond <- apply(Cret, 1:2, sum, na.rm = T) < 1  # this is a fix for low sample sizes. If Cret is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
  # fixind <- as.matrix(cbind(expand.grid(1:nsim, 1:nyears), rep(floor(maxage/3), nyears)))  # more fix
  # Cret[fixind[cond, ]] <- 1  # puts a catch in the most vulnerable age class

  # a multinomial observation model for catch-at-age data
  for (i in 1:nsim) {
    for (j in 1:nyears) {
      if (!sum( Cret[i, j,])) {
        CAA[i, j, ] <- 0
      } else {
        CAA[i, j, ] <- ceiling(-0.5 + rmultinom(1, CAA_ESS[i], Cret[i, j,]) * CAA_nsamp[i]/CAA_ESS[i])
      }
    }
  }

  # --- Simulate observed catch-at-length ----
  # a multinomial observation model for catch-at-length data
  # assumed normally-distributed length-at-age truncated at 2 standard deviations from the mean
  CAL <- array(NA, dim=c(nsim,  nyears, nCALbins))
  LFC <- rep(NA, nsim)
  vn <- (apply(N[,,,], c(1,2,3), sum) * retA[,,1:nyears]) # numbers at age that would be retained
  vn <- aperm(vn, c(1,3, 2))

  # Generate size comp data with variability in age
  tempSize <- lapply(1:nsim, genSizeCompWrap, vn, CAL_binsmid, retL, CAL_ESS, CAL_nsamp,
                     Linfarray, Karray, t0array, LenCV, truncSD=2.5)

  CAL <- aperm(array(as.numeric(unlist(tempSize, use.names=FALSE)), dim=c(nyears, length(CAL_binsmid), nsim)), c(3,1,2))


  # calculate LFC
  LFC <- unlist(lapply(tempSize, function(x) getfifth(x[nyears, ], CAL_binsmid)))
  LFC[is.na(LFC)] <- 1
  LFC[LFC<1] <- 1

  # --- Simulate index of abundance from total biomass ----
  Ierr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Isd, nyears + proyears)),
                       sdconv(1, rep(Isd, nyears + proyears))), c(nsim, nyears + proyears))
  II <- (apply(Biomass, c(1, 3), sum)^betas) * Ierr[, 1:nyears]  # apply hyperstability / hyperdepletion
  II <- II/apply(II, 1, mean)  # normalize

  # --- Calculate vulnerable and spawning biomass abundance ----
  if (nsim > 1) A <- apply(VBiomass[, , nyears, ], 1, sum)  + apply(CB[, , nyears, ], 1, sum) # Abundance before fishing
  if (nsim == 1) A <- sum(VBiomass[, , nyears, ]) +  sum(CB[,,nyears,]) # Abundance before fishing
  if (nsim > 1) Asp <- apply(SSB[, , nyears, ], 1, sum)  # SSB Abundance
  if (nsim == 1) Asp <- sum(SSB[, , nyears, ])  # SSB Abundance

  OFLreal <- A * FMSY  # the true simulated Over Fishing Limit

  # --- Simulate observed values in reference SBMSY/SB0 ----
  I3 <- apply(Biomass, c(1, 3), sum)^betas  # apply hyperstability / hyperdepletion
  I3 <- I3/apply(I3, 1, mean)  # normalize index to mean 1
  # Iref <- apply(I3[, 1:5], 1, mean) * BMSY_B0  # return the real target abundance index corresponding to BMSY
  if (nsim > 1) Iref <- apply(I3[, 1:5], 1, mean) * SSBMSY_SSB0  # return the real target abundance index corresponding to BMSY
  if (nsim == 1) Iref <- mean(I3[1:5]) * SSBMSY_SSB0

  # --- Simulate observed values in steepness ----
  if (!is.null(OM@cpars[['hsim']])) {
    hsim <- SampCpars$hsim
    hbias <- hsim/hs  # back calculate the simulated bias
    if (OM@hbiascv == 0) hbias <- rep(1, nsim)
    ObsPars$hbias <- hbias

  } else {
    if (is.null(OM@cpars[['l_hbias']])) {
      hsim <- rep(NA, nsim)
      cond <- hs > 0.6
      hsim[cond] <- 0.2 + rbeta(sum(hs > 0.6), alphaconv((hs[cond] - 0.2)/0.8, (1 - (hs[cond] - 0.2)/0.8) * OM@hbiascv),
                                betaconv((hs[cond] - 0.2)/0.8,  (1 - (hs[cond] - 0.2)/0.8) * OM@hbiascv)) * 0.8
      hsim[!cond] <- 0.2 + rbeta(sum(hs <= 0.6), alphaconv((hs[!cond] - 0.2)/0.8,  (hs[!cond] - 0.2)/0.8 * OM@hbiascv),
                                 betaconv((hs[!cond] - 0.2)/0.8, (hs[!cond] - 0.2)/0.8 * OM@hbiascv)) * 0.8
      hbias <- hsim/hs  # back calculate the simulated bias
      if (OM@hbiascv == 0) hbias <- rep(1, nsim)
      ObsPars$hbias <- hbias
    } else {
      l_hbias <- sample(OM@cpars$l_hbias, OM@nsim, replace=TRUE)

      P <- (hs-.2)/0.8
      hs_logit <- log(P/(1-P))
      P2 <- (exp(hs_logit)* l_hbias)/(1+exp(hs_logit) * l_hbias)
      P2[is.nan(P2)] <- 1
      hsim <- (P2 * 0.8) + 0.2
      hbias <- hsim/hs  # back calculate the simulated bias

      ObsPars$hbias <- hbias

    }
  }




  # --- Simulate error in observed recruitment index ----
  Recerr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(Recsd, (nyears + proyears))),
                         sdconv(1, rep(Recsd, nyears + proyears))), c(nsim, nyears + proyears))


  # --- Simulate observation error in BMSY/B0 ----
  ntest <- 20  # number of trials
  BMSY_B0bias <- array(rlnorm(nsim * ntest, mconv(1, OM@BMSY_B0biascv), sdconv(1, OM@BMSY_B0biascv)), dim = c(nsim, ntest))  # trial samples of BMSY relative to unfished
  # test <- array(BMSY_B0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0
  test <- array(SSBMSY_SSB0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0
  indy <- array(rep(1:ntest, each = nsim), c(nsim, ntest))  # index

  # indy[test > 0.9] <- NA  # interval censor
  indy[test > max(0.9, max(SSBMSY_SSB0))] <- NA  # interval censor

  BMSY_B0bias <- BMSY_B0bias[cbind(1:nsim, apply(indy, 1, min, na.rm = T))]  # sample such that BMSY_B0<90%
  ObsPars$BMSY_B0bias <- BMSY_B0bias

  # --- Implementation error time series ----

  TAC_f <- array(rlnorm(proyears * nsim, mconv(TACFrac, TACSD),
                        sdconv(TACFrac, TACSD)), c(nsim, proyears))  # composite of TAC fraction and error

  E_f <- array(rlnorm(proyears * nsim, mconv(TAEFrac, TAESD),
                      sdconv(TAEFrac, TAESD)), c(nsim, proyears))  # composite of TAC fraction and error

  SizeLim_f<-array(rlnorm(proyears * nsim, mconv(SizeLimFrac, SizeLimSD),
                          sdconv(SizeLimFrac, SizeLimSD)), c(nsim, proyears))  # composite of TAC fraction and error


  # --- Populate Data object with Historical Data ----
  Data <- new("Data", stock = "MSE")  # create a blank DLM data object
  if (reps == 1) Data <- OneRep(Data)  # make stochastic variables certain for only one rep
  Data <- replic8(Data, nsim)  # make nsim sized slots in the DLM data object
  Data@Name <- OM@Name
  Data@Year <- 1:nyears
  Data@Cat <- Cobs
  Data@Ind <- II
  Data@Rec <- apply(N[, 1, , ], c(1, 2), sum) * Recerr[, 1:nyears]
  Data@t <- rep(nyears, nsim)
  Data@AvC <- apply(Cobs, 1, mean)
  Data@Dt <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
  Data@Mort <- M * Mbias
  Data@FMSY_M <- FMSY_M * FMSY_Mbias
  # Data@BMSY_B0 <- BMSY_B0 * BMSY_B0bias
  Data@BMSY_B0 <- SSBMSY_SSB0 * BMSY_B0bias
  Data@Cref <- MSY * Crefbias
  Data@Bref <- VBMSY * Brefbias
  Data@Iref <- Iref * Irefbias
  Data@LFC <- LFC * LFCbias
  Data@LFS <- LFS[nyears,] * LFSbias
  Data@CAA <- CAA
  Data@Dep <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
  Data@Abun <- A * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
  Data@SpAbun <- Asp * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
  Data@vbK <- K * Kbias
  Data@vbt0 <- t0 * t0bias
  Data@LenCV <- LenCV # * LenCVbias
  Data@vbLinf <- Linf * Linfbias
  Data@L50 <- L50 * lenMbias
  Data@L95 <- L95 * lenMbias
  Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
  Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
  Data@steep <- hs * hbias
  Data@sigmaR <- procsd # sigmaR assumed no obs error
  Data@CAL_bins <- CAL_bins
  Data@CAL <- CAL
  MLbin <- (CAL_bins[1:(length(CAL_bins) - 1)] + CAL_bins[2:length(CAL_bins)])/2
  temp <- CAL * rep(MLbin, each = nsim * nyears)
  Data@ML <- apply(temp, 1:2, sum)/apply(CAL, 1:2, sum)
  Data@Lc <- array(MLbin[apply(CAL, 1:2, which.max)], dim = c(nsim, nyears))
  nuCAL <- CAL
  for (i in 1:nsim) for (j in 1:nyears) nuCAL[i, j, 1:match(max(1, Data@Lc[i, j]), MLbin, nomatch=1)] <- NA
  temp <- nuCAL * rep(MLbin, each = nsim * nyears)
  Data@Lbar <- apply(temp, 1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE)
  Data@MaxAge <- maxage
  Data@Units <- "unitless"
  Data@Ref <- OFLreal
  Data@Ref_type <- "Simulated OFL"
  Data@wla <- rep(a, nsim)
  Data@wlb <- rep(b, nsim)
  Data@nareas <- nareas

  # put all the operating model parameters in one table
  OMtable <- data.frame(RefY, M, Depletion, A, SSBMSY_SSB0, FMSY_M, Mgrad, Msd, procsd, Esd, dFfinal,
                        MSY=MSY, qinc, qcv, FMSY=FMSY, Linf, K, t0, hs, Linfgrad, Kgrad, Linfsd, Ksd,
                        ageM=ageM[,nyears], L5=L5[nyears, ], LFS=LFS[nyears, ], Vmaxlen=Vmaxlen[nyears, ], LFC, OFLreal,
                        Spat_targ, Size_area_1, Frac_area_1, Prob_staying, AC, L50, L95, B0, N0, SSB0, BMSY_B0,
                        TACSD,TACFrac,TAESD,TAEFrac,SizeLimSD,SizeLimFrac,Blow,MGT,
                        BMSY, SSBMSY=SSBMSY, Mexp, Fdisc,
                        LR5=LR5[nyears,], LFR=LFR[nyears,], Rmaxlen=Rmaxlen[nyears,], DR=DR[nyears,])
  OMtable <- OMtable[,order(names(OMtable))]

  Data@OM <- OMtable

  ObsTable <- as.data.frame(ObsPars)
  ObsTable <- ObsTable[,order(names(ObsTable))]
  Data@Obs <- ObsTable # put all the observation error model parameters in one table

  Data@LHYear <- OM@nyears  # Last historical year is nyears (for fixed MPs)
  histCatches <- apply(CBret, c(1, 3), sum)
  Data@MPrec <- histCatches[, nyears]
  Data@MPeff <- rep(1, nsim)
  Data@Misc <- vector("list", nsim)

  # --- Return Historical Simulations and Data from last historical year ----
  if (Hist) { # Stop the model after historical simulations are complete
    if(!silent) message("Returning historical simulations")

    nout <- t(apply(N, c(1, 3), sum))
    vb <- t(apply(VBiomass, c(1, 3), sum))
    b <- t(apply(Biomass, c(1, 3), sum))
    ssb <- t(apply(SSB, c(1, 3), sum))
    Cc <- t(apply(CB, c(1,3), sum))
    rec <- t(apply((N)[, 1, , ], c(1,2), sum))

    TSdata <- list(VB=vb, SSB=ssb, Bio=b, Catch=Cc, Rec=rec, N=nout, E_f=E_f,TAC_f=TAC_f,SizeLim_f=SizeLim_f)
    AtAge <- list(Len_age=Len_age, Wt_age=Wt_age, Sl_age=V, Mat_age=Mat_age,
                  Nage=apply(N, c(1:3), sum), SSBage=apply(SSB, c(1:3), sum), M_ageArray=M_ageArray,
                  Z=Z, FM=FM, FMret=FMret)
    MSYs <- list(MSY=MSY, FMSY=FMSY, VBMSY=VBMSY, UMSY=UMSY,
                 SSBMSY=SSBMSY, BMSY_B0=BMSY_B0, SSBMSY_SSB0=SSBMSY_SSB0,
                 SSB0=SSB0, B0=B0)

    StockPars$Depletion <- Depletion
    FleetPars$qs <- qs
    SampPars <- c(StockPars, FleetPars, ObsPars, ImpPars)
    Data@Misc <- list()
    HistData <- list(SampPars=SampPars, TSdata=TSdata, AtAge=AtAge, MSYs=MSYs,
                     Data=Data)
    return(HistData)
  }

  # assign('Data',Data,envir=.GlobalEnv) # for debugging fun

  # --- Check MPs ----
  if (is.na(MPs[1])) CheckMPs <- TRUE
  if (CheckMPs) {
    if(!silent) message("Determining available methods")  # print an progress report
    PosMPs <- Can(Data, timelimit = timelimit)  # list all the methods that could be applied to a DLM data object
    if (is.na(MPs[1])) {
      MPs <- PosMPs  # if the user does not supply an argument MPs run the MSE for all available methods
      if(!silent) message("No MPs specified: running all available")
    }
    if (!is.na(MPs[1])) {
      cant <- MPs[!MPs %in% PosMPs]
      if (length(cant) > 0) {
        if(!silent) message("Cannot run some MPs: ")
        if(!silent) print(DLMdiag(Data, "not available", funcs1=cant, timelimit = timelimit))
      }
      MPs <- MPs[MPs %in% PosMPs]  # otherwise run the MSE for all methods that are deemed possible
    }
    if (length(MPs) == 0) {
      if(!silent) message(Cant(Data, timelimit = timelimit))
      stop("MSE stopped: no viable methods \n\n")  # if none of the user specied methods are possible stop the run
    }
  }

  ok <- rep(TRUE, length(MPs))
  for (mm in seq_along(MPs)) {
    test <- try(get(MPs[mm]), silent=TRUE)
    if (!class(test) == 'MP') {
      ok[mm] <- FALSE
      if (class(test) == 'try-error') {
        message('Object ', paste(MPs[mm], ""), " does not exist - Ignoring")
      } else message('Dropping MP: ', paste(MPs[mm], ""), " - Not class 'MP'")
    }

  }

  MPs <- MPs[ok]

  nMP <- length(MPs)  # the total number of methods used

  if (nMP < 1) stop("No valid MPs found")
  MSElist <- list(Data)[rep(1, nMP)]  # create a data object for each method (they have identical historical data and branch in projected years)

  B_BMSYa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected B_BMSY
  F_FMSYa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected F_FMSY
  Ba <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected Biomass
  SSBa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected SSB
  VBa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected vulnerable biomass
  FMa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected fishing mortality rate
  Ca <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected removed catch
  CaRet <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected retained catch
  TACa <- array(NA, dim = c(nsim, nMP, proyears))  # store the projected TAC recommendation
  Effort <- array(NA, dim = c(nsim, nMP, proyears))  # store the Effort
  PAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the population-at-age in last projection year
  CAAout <- array(NA, dim = c(nsim, nMP, maxage))  # store the catch-at-age in last projection year
  CALout <- array(NA, dim = c(nsim, nMP, nCALbins))  # store the population-at-length in last projection year

  # SPRa <- array(NA,dim=c(nsim,nMP,proyears)) # store the Spawning Potential Ratio

  # --- Calculate MSY statistics for each projection year ----
  MSY_P <- array(MSY, dim=c(nsim, nMP, proyears))
  FMSY_P <- array(FMSY, dim=c(nsim, nMP, proyears))
  SSBMSY_P <- array(SSBMSY, dim=c(nsim, nMP, proyears))

  if (AnnualMSY) {

    if(!silent) message("Calculating MSY reference points for each projection year")
    for (y in 1:proyears) {
      if(!silent) cat('.')
      if (!silent) flush.console()
      MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age, V, maxage, R0, SRrel, hs, yr=nyears+y)
      MSY_P[,,y] <- MSYrefsYr[1, ]
      FMSY_P[,,y] <- MSYrefsYr[2,]
      SSBMSY_P[,,y] <- MSYrefsYr[3,]

    }
    if(!silent) cat("\n")
  }

  if (length(interval) != nMP) interval <- rep(interval, nMP)[1:nMP]
  if (!all(interval == interval[1])) {
    message("Variable management intervals:")
    df <- data.frame(MP=MPs,interval=interval)
    message(paste(capture.output(print(df)), collapse = "\n"))
  }

  # --- Begin loop over MPs ----
  mm <- 1 # for debugging

  for (mm in 1:nMP) {  # MSE Loop over methods

    if(!silent) message(mm, "/", nMP, " Running MSE for ", MPs[mm])  # print a progress report
    checkNA <- NA # save number of NAs

    # reset selectivity parameters for projections
    L5_P <- L5
    LFS_P <- LFS
    Vmaxlen_P <- Vmaxlen
    SLarray_P <- SLarray # selectivity at length array - projections
    V_P <- V  #  selectivity at age array - projections

    # reset retention parametersfor projections
    LR5_P <- LR5
    LFR_P <- LFR
    Rmaxlen_P <- Rmaxlen
    retA_P <- retA # retention at age array - projections
    retL_P <- retL # retention at length array - projections

    Fdisc_P <- Fdisc # Discard mortality for projectons
    DR_P <- DR # Discard ratio for projections

    # projection arrays
    N_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    Biomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    VBiomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSN_P <-array(NA, dim = c(nsim, maxage, proyears, nareas))
    SSB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    FM_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    FM_Pret <- array(NA, dim = c(nsim, maxage, proyears, nareas)) # retained F
    FM_nospace <- array(NA, dim = c(nsim, maxage, proyears, nareas))  # stores prospective F before reallocation to new areas
    FML <- array(NA, dim = c(nsim, nareas))  # last apical F
    Z_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    CB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
    CB_Pret <- array(NA, dim = c(nsim, maxage, proyears, nareas)) # retained catch

    # indexes
    SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
    SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, 1 + nyears, 1:nareas))  # Trajectory year
    SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, 1, 1:nareas))
    SYt <- SAYRt[, c(1, 3)]
    SAYt <- SAYRt[, 1:3]
    SR <- SAYR[, c(1, 4)]
    SA1 <- SAYR[, 1:2]
    S1 <- SAYR[, 1]
    SY1 <- SAYR[, c(1, 3)]
    SAY1 <- SAYRt[, 1:3]
    SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
    SY <- SYA[, 1:2]
    SA <- SYA[, c(1, 3)]
    SAY <- SYA[, c(1, 3, 2)]
    S <- SYA[, 1]

    # -- First projection year ----
    y <- 1
    NextYrN <- lapply(1:nsim, function(x)
      popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB[x,,nyears, ]), Ncurr=N[x,,nyears,],
                     Zcurr=Z[x,,nyears,], PerrYr=Perr_y[x, nyears+maxage-1], hs=hs[x],
                     R0a=R0a[x,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                     mov=mov[x,,,], SRrel=SRrel[x]))

    N_P[,,1,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2))

    Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAY1]  # Calculate biomass
    VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
    SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAY1]  # Calculate spawning stock numbers
    SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAY1]
    FML <- apply(FM[, , nyears, ], c(1, 3), max)

    # -- apply MP in initial projection year ----
    # Combined MP ----
    runMP <- applyMP(MSElist[[mm]], MPs = MPs[mm], reps = reps, silent=TRUE)  # Apply MP

    MPRecs <- runMP[[1]][[1]] # MP recommendations
    Data_p <- runMP[[2]] # Data object object with saved info from MP
    Data_p@TAC <- MPRecs$TAC

    # calculate pstar quantile of TAC recommendation dist
    TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, na.rm = T)
    checkNA[y] <- sum(is.na(TACused))
    LastEi <- rep(1,nsim) # no effort adjustment
    LastSpatial <- array(MPA[nyears,], dim=c(nareas, nsim)) #
    LastAllocat <- rep(1, nsim) # default assumption of reallocation of effort to open areas
    LastCatch <- apply(CB[,,nyears,], 1, sum)

    MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, nsim,
                              LastEi, LastSpatial, LastAllocat, LastCatch,
                              TACused, maxF,
                              LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                              L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                              Fdisc_P, DR_P,
                              M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                              TAC_f, E_f, SizeLim_f,
                              VBiomass_P, Biomass_P, FinF, Spat_targ,
                              CAL_binsmid, Linf, Len_age, maxage, nareas, Asize, nCALbins,
                              qs, qvar, qinc)

    TACa[, mm, y] <- MPCalcs$TACrec # recommended TAC
    LastSpatial <- MPCalcs$Si
    LastAllocat <- MPCalcs$Ai
    LastEi <- MPCalcs$Ei # adjustment to effort
    LastCatch <- MPCalcs$TACrec

    Effort[, mm, y] <- MPCalcs$Effort
    CB_P <- MPCalcs$CB_P # removals
    CB_Pret <- MPCalcs$CB_Pret # retained catch
    FM_P <- MPCalcs$FM_P # fishing mortality
    FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality
    Z_P <- MPCalcs$Z_P # total mortality

    #### DEBUG ####

    # fs <- -log(1 - apply(CB_P[,,y,], c(1), sum, na.rm=TRUE)/apply(VBiomass_P[,,y,]+CB_P[,,y,], c(1), sum, na.rm=TRUE))
    # fs/FMSY_P[,mm,y]

    #### DEBUG ####

    retA_P <- MPCalcs$retA_P # retained-at-age

    retL_P <- MPCalcs$retL_P # retained-at-length
    V_P <- MPCalcs$V_P  # vulnerable-at-age
    SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length

    upyrs <- 1 + (0:(floor(proyears/interval[mm]) - 1)) * interval[mm]  # the years in which there are updates (every three years)
    if(!silent) {
      cat(".")
      flush.console()
    }

    # --- Begin projection years ----
    for (y in 2:proyears) {
      if(!silent) {
        cat(".")
        flush.console()
      }

      SelectChanged <- FALSE
      if (AnnualMSY) {
        if (any(range(retA_P[,,nyears+y] - retA[,,nyears+y]) !=0)) SelectChanged <- TRUE
        if (any(range(V_P[,,nyears+y] - V[,,nyears+y]) !=0))  SelectChanged <- TRUE
      }

      # -- Calculate MSY stats for this year ----
      if (AnnualMSY & SelectChanged) { #
        MSYrefsYr <- sapply(1:nsim, optMSY_eq, M_ageArray, Wt_age, Mat_age,
                            V_P, maxage, R0, SRrel, hs, yr=nyears+y)
        MSY_P[,mm,y] <- MSYrefsYr[1, ]
        FMSY_P[,mm,y] <- MSYrefsYr[2,]
        SSBMSY_P[,mm,y] <- MSYrefsYr[3,]

      }

      TACa[, mm, y] <- TACa[, mm, y-1] # TAC same as last year unless changed

      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
      SAYt <- SAYRt[, 1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt <- SAYRt[, c(1, 3)]
      SAY1R <- as.matrix(expand.grid(1:nsim, 1:maxage, y - 1, 1:nareas))
      SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
      SY <- SAYR[, c(1, 3)]
      SA <- SAYR[, 1:2]
      S1 <- SAYR[, 1]

      SAY <- SAYR[, 1:3]
      S <- SAYR[, 1]
      SR <- SAYR[, c(1, 4)]
      SA2YR <- as.matrix(expand.grid(1:nsim, 2:maxage, y, 1:nareas))
      SA1YR <- as.matrix(expand.grid(1:nsim, 1:(maxage - 1), y -1, 1:nareas))

      NextYrN <- lapply(1:nsim, function(x)
        popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB_P[x,,y-1, ]), Ncurr=N_P[x,,y-1,],
                       Zcurr=Z_P[x,,y-1,], PerrYr=Perr_y[x, y+nyears+maxage-1], hs=hs[x],
                       R0a=R0a[x,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                       mov=mov[x,,,], SRrel=SRrel[x]))

      N_P[,,y,] <- aperm(array(unlist(NextYrN), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2))
      Biomass_P[SAYR] <- N_P[SAYR] * Wt_age[SAYt]  # Calculate biomass
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]  # Calculate vulnerable biomass
      SSN_P[SAYR] <- N_P[SAYR] * Mat_age[SAYt]  # Calculate spawning stock numbers
      SSB_P[SAYR] <- SSN_P[SAYR] * Wt_age[SAYt]  # Calculate spawning stock biomass

      # --- An update year ----
      if (y %in% upyrs) {
        # rewrite the DLM object and run the TAC function
        yind <- upyrs[match(y, upyrs) - 1]:(upyrs[match(y, upyrs)] - 1)

        # use the retained catch
        CBtemp <-  CB_Pret[, , yind, , drop=FALSE] # retained catch-at-age
        CNtemp <- retA_P[,,yind+nyears, drop=FALSE] * apply(N_P[,,yind,, drop=FALSE], c(1,2,3), sum) # retained age structure

        CBtemp[is.na(CBtemp)] <- tiny
        CBtemp[!is.finite(CBtemp)] <- tiny
        CNtemp[is.na(CNtemp)] <- tiny
        CNtemp[!is.finite(CNtemp)] <- tiny

        Cobs <- Cbiasa[, nyears + yind] * Cerr[, nyears + yind] * apply(CBtemp, c(1, 3), sum, na.rm = T)
        Cobs[is.na(Cobs)] <- tiny
        Recobs <- Recerr[, nyears + yind] * apply(array(N_P[, 1, yind, ], c(nsim, interval[mm], nareas)), c(1, 2), sum)

        CAA <- array(0, dim = c(nsim, interval[mm], maxage))  # Catch  at age array

        # # a multinomial observation model for catch-at-age data
        for (i in 1:nsim) {
          for (j in 1:interval[mm]) {
            # if (all(CNtemp[i, , j]<1)) { # this is a fix for low sample sizes. If CAA is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
            #   CNtemp[i, floor(maxage/3), j] <- 1
            # }
            CAA[i, j, ] <- ceiling(-0.5 + rmultinom(1, CAA_ESS[i], CNtemp[i, , j]) * CAA_nsamp[i]/CAA_ESS[i])   # a multinomial observation model for catch-at-age data
          }
        }


        ## Calculate CAL ####
        CAL <- array(NA, dim = c(nsim, interval[mm], nCALbins))  # the catch at length array
        # # a multinomial observation model for catch-at-length data

        vn <- (apply(N_P[,,,], c(1,2,3), sum) * retA_P[,,(nyears+1):(nyears+proyears)]) # numbers at age that would be retained
        vn <- aperm(vn, c(1,3,2))

        nyrs <- length(yind)
        tempSize <- lapply(1:nsim, genSizeCompWrap, vn[,yind,, drop=FALSE], CAL_binsmid, retL_P,
                           CAL_ESS, CAL_nsamp,
                           Linfarray[,nyears + yind, drop=FALSE],
                           Karray[,nyears + yind, drop=FALSE],
                           t0array[,nyears + yind,drop=FALSE], LenCV, truncSD=2)
        CAL <- aperm(array(as.numeric(unlist(tempSize, use.names=FALSE)), dim=c(length(yind), length(CAL_binsmid), nsim)), c(3,1,2))

        # check CAL
        # chkCAL <- apply(CAL, c(1,2), sum)
        # if (any(chkCAL==0)) {
        #   for (ss in 1:nrow(chkCAL)) {
        #     for (yy in 1:ncol(chkCAL)) {
        #       if (chkCAL[ss,yy]==0) {
        #         bn <- ceiling(c(0.5*nCALbins, 0.5*nCALbins+1))
        #         CAL[ss,yy,bn] <- 1 # if CAL is empty - put 1 individual in two middle bins - for robustness
        #       }
        #     }
        #   }
        # }

        # calculate LFC - approx 5th percentile
        LFC <- unlist(lapply(tempSize, function(x) getfifth(x[nrow(x), ], CAL_binsmid)))
        LFC[is.na(LFC)] <- 1
        LFC[LFC<1] <- 1

        I2 <- (cbind(apply(Biomass, c(1, 3), sum), apply(Biomass_P, c(1, 3), sum)[, 1:(y - 1)])^betas) * Ierr[, 1:(nyears + (y - 1))]
        I2[is.na(I2)] <- tiny
        I2 <- I2/apply(I2, 1, mean)

        Depletion <- apply(SSB_P[, , y, ], 1, sum)/SSB0 # apply(SSB[, , 1, ], 1, sum)
        Depletion[Depletion < tiny] <- tiny

        NextYrNtemp <- lapply(1:nsim, function(x)
          popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB_P[x,,y, ]), Ncurr=N_P[x,,y,],
                         Zcurr=matrix(M_ageArray[x,,y+nyears], ncol=nareas, nrow=maxage),
                         PerrYr=Perr_y[x, y+nyears+maxage-1], hs=hs[x],
                         R0a=R0a[x,], SSBpR=SSBpR[x,], aR=aR[x,], bR=bR[x,],
                         mov=mov[x,,,], SRrel=SRrel[x]))

        N_PNext <- aperm(array(unlist(NextYrNtemp), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2))
        VBiomassNext <- VBiomass_P
        VBiomassNext[SAYR] <- N_PNext * Wt_age[SAYt] * V_P[SAYt]  # Calculate vulnerable for abundance

        A <- apply(VBiomassNext[, , y, ], 1, sum)
        # A <- apply(VBiomass_P[, , y, ], 1, sum)

        A[is.na(A)] <- tiny
        Asp <- apply(SSB_P[, , y, ], 1, sum)  # SSB Abundance
        Asp[is.na(Asp)] <- tiny
        OFLreal <- A * FMSY_P[,mm,y]

        # - update data object ----
        # assign all the new data
        MSElist[[mm]]@OM$A <- A
        MSElist[[mm]]@Year <- 1:(nyears + y - 1)
        MSElist[[mm]]@Cat <- cbind(MSElist[[mm]]@Cat, Cobs)
        MSElist[[mm]]@Ind <- I2
        MSElist[[mm]]@Rec <- cbind(MSElist[[mm]]@Rec, Recobs)
        MSElist[[mm]]@t <- rep(nyears + y, nsim)
        MSElist[[mm]]@AvC <- apply(MSElist[[mm]]@Cat, 1, mean)
        MSElist[[mm]]@Dt <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
        oldCAA <- MSElist[[mm]]@CAA
        MSElist[[mm]]@CAA <- array(0, dim = c(nsim, nyears + y - 1, maxage))
        MSElist[[mm]]@CAA[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAA[, 1:(nyears + y - interval[mm] - 1), ] # there is some bug here sometimes oldCAA (MSElist[[mm]]@CAA previously) has too many years of observations
        MSElist[[mm]]@CAA[, nyears + yind, ] <- CAA
        MSElist[[mm]]@Dep <- Dbias * Depletion * rlnorm(nsim, mconv(1, Derr), sdconv(1, Derr))
        MSElist[[mm]]@Abun <- A * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
        MSElist[[mm]]@SpAbun <- Asp * Abias * rlnorm(nsim, mconv(1, Aerr), sdconv(1, Aerr))
        MSElist[[mm]]@CAL_bins <- CAL_bins
        oldCAL <- MSElist[[mm]]@CAL
        MSElist[[mm]]@CAL <- array(0, dim = c(nsim, nyears + y - 1, nCALbins))
        MSElist[[mm]]@CAL[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAL[, 1:(nyears + y - interval[mm] - 1), ]# there is some bug here: sometimes oldCAL (MSElist[[mm]]@CAL previously) has too many years of observations
        MSElist[[mm]]@CAL[, nyears + yind, ] <- CAL[, 1:interval[mm], ]

        temp <- CAL * rep(MLbin, each = nsim * interval[mm])
        MSElist[[mm]]@ML <- cbind(MSElist[[mm]]@ML, apply(temp, 1:2, sum)/apply(CAL, 1:2, sum))
        MSElist[[mm]]@Lc <- cbind(MSElist[[mm]]@Lc, array(MLbin[apply(CAL, 1:2, which.max)], dim = c(nsim, interval[mm])))
        nuCAL <- CAL
        for (i in 1:nsim) for (j in 1:interval[mm]) nuCAL[i, j, 1:match(max(1, MSElist[[mm]]@Lc[i, j]), MLbin,nomatch=1)] <- NA
        temp <- nuCAL * rep(MLbin, each = nsim * interval[mm])
        MSElist[[mm]]@Lbar <- cbind(MSElist[[mm]]@Lbar, apply(temp,1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE))

        MSElist[[mm]]@LFC <- LFC * LFCbias
        MSElist[[mm]]@LFS <- LFS[nyears + y,] * LFSbias

        # update growth, maturity estimates for current year
        MSElist[[mm]]@vbK <-  Karray[, nyears+y] * Kbias
        MSElist[[mm]]@vbt0 <- t0 * t0bias

        MSElist[[mm]]@vbLinf <- Linfarray[, nyears+y] * Linfbias
        MSElist[[mm]]@L50 <- L50array[, nyears+y] * lenMbias
        MSElist[[mm]]@L95 <- L95array[, nyears+y] * lenMbias
        MSElist[[mm]]@L95[is.na(MSElist[[mm]]@L95)]<-MSElist[[mm]]@vbLinf # this is just to robustify 'numbers models' like Grey Seal that do not generate (and will never use) real length observations
        MSElist[[mm]]@L95[MSElist[[mm]]@L95 > 0.9 * MSElist[[mm]]@vbLinf] <- 0.9 * MSElist[[mm]]@vbLinf[MSElist[[mm]]@L95 > 0.9 * MSElist[[mm]]@vbLinf]  # Set a hard limit on ratio of L95 to Linf
        MSElist[[mm]]@L50[MSElist[[mm]]@L50 > 0.9 * MSElist[[mm]]@L95] <- 0.9 * MSElist[[mm]]@L95[MSElist[[mm]]@L50 > 0.9 * MSElist[[mm]]@L95]  # Set a hard limit on ratio of L95 to Linf

        MSElist[[mm]]@Ref <- OFLreal
        MSElist[[mm]]@Ref_type <- "Simulated OFL"
        MSElist[[mm]]@Misc <- Data_p@Misc
        MSElist[[mm]]@MPrec <- MPCalcs$TACrec # last MP  TAC recommendation
        MSElist[[mm]]@MPeff <- Effort[, mm, y-1] # last recommended effort

        # assign('Data',MSElist[[mm]],envir=.GlobalEnv) # for debugging fun

        # apply combined MP ----
        if("DataOut"%in%names(control))if(control$DataOut == y) return(MSElist)

        runMP <- applyMP(MSElist[[mm]], MPs = MPs[mm], reps = reps, silent=TRUE)  # Apply MP

        MPRecs <- runMP[[1]][[1]] # MP recommendations
        Data_p <- runMP[[2]] # Data object object with saved info from MP
        Data_p@TAC <- MPRecs$TAC

        # calculate pstar quantile of TAC recommendation dist
        TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, na.rm = T)

        MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, nsim,
                                  LastEi, LastSpatial, LastAllocat, LastCatch,
                                  TACused, maxF,
                                  LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                  L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                  Fdisc_P, DR_P,
                                  M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                  TAC_f, E_f, SizeLim_f,
                                  VBiomass_P, Biomass_P, FinF, Spat_targ,
                                  CAL_binsmid, Linf, Len_age, maxage, nareas, Asize,  nCALbins,
                                  qs, qvar, qinc)

        TACa[, mm, y] <- TACused # MPCalcs$TACrec # recommended TAC
        LastSpatial <- MPCalcs$Si
        LastAllocat <- MPCalcs$Ai
        LastEi <- MPCalcs$Ei # adjustment to effort

        Effort[, mm, y] <- MPCalcs$Effort


        CB_P <- MPCalcs$CB_P # removals
        CB_Pret <- MPCalcs$CB_Pret # retained catch

        LastCatch <- apply(CB_Pret[,,y,], 1, sum, na.rm=TRUE)

        FM_P <- MPCalcs$FM_P # fishing mortality
        FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality
        Z_P <- MPCalcs$Z_P # total mortality

        retA_P <- MPCalcs$retA_P # retained-at-age
        retL_P <- MPCalcs$retL_P # retained-at-length
        V_P <- MPCalcs$V_P  # vulnerable-at-age
        SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length


      } else {
        # --- Not an update yr ----

        NoMPRecs <- MPRecs
        NoMPRecs[lapply(NoMPRecs, length) > 0 ] <- NULL
        NoMPRecs$Spatial <- NA

        MPCalcs <- CalcMPDynamics(NoMPRecs, y, nyears, proyears, nsim,
                                  LastEi, LastSpatial, LastAllocat, LastCatch,
                                  TACused, maxF,
                                  LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P,
                                  L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P,
                                  Fdisc_P, DR_P,
                                  M_ageArray, FM_P, FM_Pret, Z_P, CB_P, CB_Pret,
                                  TAC_f, E_f, SizeLim_f,
                                  VBiomass_P, Biomass_P, FinF, Spat_targ,
                                  CAL_binsmid, Linf, Len_age, maxage, nareas, Asize,  nCALbins,
                                  qs, qvar, qinc)

        TACa[, mm, y] <- TACused #  MPCalcs$TACrec # recommended TAC
        LastSpatial <- MPCalcs$Si
        LastAllocat <- MPCalcs$Ai
        LastEi <- MPCalcs$Ei
        Effort[, mm, y] <- MPCalcs$Effort
        CB_P <- MPCalcs$CB_P # removals
        CB_Pret <- MPCalcs$CB_Pret # retained catch

        LastCatch <- apply(CB_Pret[,,y,], 1, sum, na.rm=TRUE)

        FM_P <- MPCalcs$FM_P # fishing mortality
        FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality
        Z_P <- MPCalcs$Z_P # total mortality

        retA_P <- MPCalcs$retA_P # retained-at-age
        retL_P <- MPCalcs$retL_P # retained-at-length
        V_P <- MPCalcs$V_P  # vulnerable-at-age
        SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length

      }  # not an update year
      checkNA[y] <- sum(is.na(TACused))
    }  # end of year

    B_BMSYa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE)/SSBMSY_P[,mm,]  # SSB relative to SSBMSY

    FMa[, mm, ] <- -log(1 - apply(CB_P, c(1, 3), sum, na.rm=TRUE)/apply(VBiomass_P+CB_P, c(1, 3), sum, na.rm=TRUE))
    F_FMSYa[, mm, ] <- FMa[, mm, ]/FMSY_P[,mm,]

    Ba[, mm, ] <- apply(Biomass_P, c(1, 3), sum, na.rm=TRUE) # biomass
    SSBa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm=TRUE) # spawning stock biomass
    VBa[, mm, ] <- apply(VBiomass_P, c(1, 3), sum, na.rm=TRUE) # vulnerable biomass

    Ca[, mm, ] <- apply(CB_P, c(1, 3), sum, na.rm=TRUE) # removed
    CaRet[, mm, ] <- apply(CB_Pret, c(1, 3), sum, na.rm=TRUE) # retained catch

    # Store Pop and Catch-at-age and at-length for last projection year
    PAAout[ , mm, ] <- apply(N_P[ , , proyears, ], c(1,2), sum) # population-at-age

    CNtemp <- apply(CB_Pret, c(1,2,3), sum)/Wt_age[(nyears+1):nyears+proyears]
    CAAout[ , mm, ] <- CNtemp[,,proyears] # nsim, maxage # catch-at-age
    CALout[ , mm, ] <- CAL[,max(dim(CAL)[2]),] # catch-at-length in last year

    if (!silent) {
      cat("\n")
      if (all(checkNA != nsim) & !all(checkNA == 0)) {
        # print number of NAs
        # message(checkNA)
        # message(checkNA[upyrs])
        ntot <- sum(checkNA[upyrs])
        totyrs <- sum(checkNA[upyrs] >0)
        nfrac <- round(ntot/(length(upyrs)*nsim),2)*100

        message(totyrs, ' years had TAC = NA for some simulations (', nfrac, "% of total simulations)")
        message('Used TAC_y = TAC_y-1')
      }

    }

    if (!parallel)
      if("progress"%in%names(control))
        if(control$progress)
          shiny::incProgress(1/nMP, detail = round(mm*100/nMP))

  }  # end of mm methods

  # Miscellaneous reporting
  if(PPD)Misc<-MSElist

  # rescale effort to today

  Effort <-  Effort/array(FinF, dim=c(nsim, nMP, proyears))

  ## Create MSE Object ####
  MSEout <- new("MSE", Name = OM@Name, nyears, proyears, nMPs=nMP, MPs, nsim,
                Data@OM, Obs=Data@Obs, B_BMSY=B_BMSYa, F_FMSY=F_FMSYa, B=Ba,
                SSB=SSBa, VB=VBa, FM=FMa, CaRet, TAC=TACa, SSB_hist = SSB, CB_hist = CB,
                FM_hist = FM, Effort = Effort, PAA=PAAout, CAA=CAAout, CAL=CALout, CALbins=CAL_binsmid,
                Misc = Misc)
  # Store MSE info
  attr(MSEout, "version") <- packageVersion("DLMtool")
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version

  MSEout
}




