
# ====================================================================
# === CASAL2OM =======================================================
# ====================================================================

# April 2021

#' Reads MLE estimates from CASAL file structure into an operating model
#'
#'
#' @description A function that uses the file location of a fitted CASAL assessment model including input files to population the
#' various slots of an operating model with MLE parameter estimates. The function mainly populates the Stock and Fleet portions
#' of the operating model; the user still needs to parameterize most of the observation and implementation portions of the operating model.
#' @param CASALdir A folder with Stock Synthesis input and output files in it
#' @param Obs The observation model (class Obs).
#' @param Imp The implementation model (class Imp).
#' @param Name The common name of the operating model
#' @param Agency The fishery management agency
#' @param Region The geographical location
#' @param Sponsor Who funded the work
#' @param Latitude In degrees north
#' @param Longitude In degrees west
#' @param nsim The number of simulations to take for parameters with uncertainty (for OM@@cpars custom parameters).
#' @param proyears The number of projection years for MSE
#' @param interval The number of years between management updates
#' @param pstar The quantile for TAC management given stochasticity
#' @param maxF The maximum allowable F in the operating model.
#' @param reps The number of stochastic replicates within each simulation in the operating model.
#' @param seed The random seed for the operating model.
#' @param Common_Name The name of the species
#' @param Species The species latin name
#' @param Source Reference to assessment documentation e.g. a url
#' @param Author Who did the assessment
#' @param ... Arguments to pass to \link[r4ss]{SS_output}.
#' @return An object of class OM.
#' @author T. Carruthers
#' @export
#' @seealso \link{SS2OM}


#' @importFrom stats acf

CASAL2OM<-function(CASALdir,Obs=DLMtool::Precise_Unbiased, Imp=DLMtool::Perfect_Imp,
                   Name=NA,Agency=NA, Region=NA, Sponsor=NA, Latitude=NA, 
                   Longitude=NA, nsim=48,proyears=50, interval=4, pstar=0.5, maxF=2, 
                   reps=1,seed=1, Source=NA,
                    
                  Common_Name=NA,Species=NA, Source=NA,Author=NA){
  
  # Known issues / areas for improvement
  # - currently MLE only
  # - the maturity ogive ripping only works for type "allvalues_bounded" - user defined minage and maxage assign a vector that long assuming all values below minage are 0 and all above maxage are 1
  # - the mortality rate is only set up here for a single-age invariant M
  # - no sex specific population dynamics are currently included
  
  
  files<-list.files(CASALdir)
  runname<-strsplit(files[grep("casal_estimation",files)],"_")[[1]][1]
  
  est <- readLines(paste0(CASALdir,"/",files[grep("casal_estimation",files)])) # estimation set up
  out <- readLines(paste0(CASALdir,"/",files[grep("casal_output",files)]))     # what outputs to return
  pop <- readLines(paste0(CASALdir,"/",files[grep("casal_population",files)])) # input data and parameters
  mpd <- readLines(paste0(CASALdir,"/",files[grep("_mpd",files)]))             # q B0 recruitment and selectivity  
  output <- readLines(paste0(CASALdir,"/",files[grep("_output.log",files)]))
  runobj <- readLines(paste0(CASALdir,"/",files[grep("objectives",files)])) # MCMC objective function output
  runout <- readLines(paste0(CASALdir,"/",files[grep(paste("output",runname,"casal_run",sep="_"),files)])) #
  runsamp <- readLines(paste0(CASALdir,"/",files[grep("samples",files)]))
  
  rip<-function(obj,lookup,pos=2,sep=" ",relline=0){
    if(pos!="all"){
      outy<-strsplit(obj[grep(lookup,obj)+relline],sep)[[1]][pos]  # get a value from obj, identified by lookup, in postion pos in the vector, separated by sep, rellines below lookup
    }else{
      outy<-strsplit(obj[grep(lookup,obj)+relline],sep)[[1]]  # get a value from obj, identified by lookup, in postion pos in the vector, separated by sep, rellines below lookup
    }
    outy
  }
  
  CASALpars<-function(est,output,out,pop){
    
    est_type <- rip(est, "@estimator")
    Author   <- rip(output,"User name",sep=": ") 
    years    <- rip(out,"years","all") 
    years    <- years[2:length(years)]
    nyears   <- length(years)
    maxage   <- as.integer(rip(pop,"@max_age"))
    SR_type  <- rip(pop,"SR")
    Growth_type <- rip(pop,"@size_at_age_type von_Bert")
    R0       <- as.numeric(rip(output,"* R0",pos=1,relline=1)) 
    Perr     <- as.numeric(rip(pop,"sigma_r"))
    h       <- as.numeric(rip(pop,"steepness"))
    
    # if all ages
    M        <- as.numeric(rip(pop,"@natural_mortality",relline=1))
    Linf     <- as.numeric(rip(pop,"Linf"))
    t0       <- as.numeric(rip(pop,"t0"))
    K        <- as.numeric(rip(pop,"k"))
    LenCV    <- as.numeric(rip(pop,"cv"))
    a        <- as.numeric(rip(pop,"a "))
    b        <- as.numeric(rip(pop,"b "))
    mat_age<-rep(1,maxage)
   
    # if @maturity_props
    matinst<- rip(pop,"@maturity_props",pos="all",relline=1)
    mat_age[matinst[3]:matinst[4]]<-matinst[5:length(matinst)]
    mat_age[1:matinst[3]]<-0
    mat_age<-as.numeric(mat_age)
    
    pars     <- list(est_type=est_type, Author=Author, years=years, 
                     nyears=nyears, maxage=maxage, SR_type=SR_type,Growth_type=Growth_type,
                     R0=R0, Perr=Perr, h=h, M=M, Linf=Linf, t0=t0, K=K, LenCV=LenCV, 
                     a=a, b=b, mat_age=mat_age)
    pars
    
  }
  
 
  # Construct OM
  Stock <- new("Stock")
  Fleet <- new("Fleet")
  OM <- new("OM", Stock = Stock, Fleet = Fleet, Obs = Obs, Imp = Imp)
  
  # OM Misc
  OM@Name <- Name
  OM@Agency <- Agency
  OM@Region <- Region
  OM@Sponsor <- Sponsor
  OM@Latitute <- Latitude
  OM@Longitude <- Longitude
  OM@nsim <- nsim
  OM@proyears <- proyears
  OM@interval <- interval
  OM@pstar <- pstar
  OM@maxF <- maxF
  OM@reps <- reps
  #cpars
  OM@seed <- seed
  OM@source <- Source
  OM@Common_Name <- Common_Name
  OM@Species <- Species
  
  # CASAL pars
  pars<-CASALpars(est,output,out,pop)
  OM@maxage<-pars$maxage
  OM@R0<-pars$R0
  OM@M<-rep(pars$M,2)
  OM@Msd<-c(0.005,0.01)
  OM@h<-rep(pars$h,2)
  OM@SRrel<-1 # match(pars$SR_type,c("BH","")) # don't know how CASAL reports the Ricker, maybe just 'R'?
  OM@Perr<-rep(pars$Perr,2)
  #OM@AC 
  OM@Linf <-rep(pars$Linf,2)
  OM@K<-rep(pars$K,2)
  OM@t0 <- rep(pars$t0, 2)
  OM@LenCV <-rep(pars$LenCV,2)
  #L50
  #L50_95
  
  # YOU GOT TO HERE 
  OM@D<-"ermm"!!!
  OM@a<-pars$a
  OM@b<-pars$b
  OM@Size_area_1<-OM@Frac_area_1<-c(0.1,0.1)
  OM@Prob_staying<-c(0.5,0.5)
  #OM@Fdist<-
  OM@nyears<-pars$nyears

  
  
  
  
  
}


CASALdir="C:/Users/tcarruth/Documents/GitHub/DLMDev/Case_Studies/Z - INCOMPLETE/Toothfish_Patagonia_CCAMLR/data/CASAL/"
nsim<-4
proyears<-50
reps=1
maxF=2
seed=1
Obs=DLMtool::Precise_Unbiased
Imp=DLMtool::Perfect_Imp
Name=NA
Source=NA
Author=NA


read






