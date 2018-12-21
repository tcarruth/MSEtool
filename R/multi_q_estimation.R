

#' optimize for catchability (q) and fishing dist for a MICE model
#'
#' Function optimizes catchability (q, where F=qE) required to get to user-specified stock
#' depletion across stocks and fleets if there are relationships among stocks
#'
#' @param x Integer, the simulation number
#' @param StockPars A list of sampled stock parameters, one list element per stock
#' @param FleetPars A hierarcical list of sampled fleet parameters, first list level is stock, second is fleet
#' @param np The number of stocks
#' @param nf The number of fleets
#' @param nareas The number of areas
#' @param maxage The maximum number of modelled ages
#' @param nyears The number of historical 'spool-up' years (from unfished to now)
#' @param N An array of stock numbers [nsim,np,maxage,nyears,nareas] - only the values from the first year are used
#' @param VF An array of vulnerability [nsim,np,nf,maxage,nyears+proyears]
#' @param FretA An array of retention [nsim,np,nf,maxage,nyears+proyears]
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param MPA An of spatial closures by year [np,nf,nyears+proyears,nareas]
#' @param CatchFrac A list of stock-specific fleet fractions of current catch list[[stock]][nsim, nf]
#' @param bounds Bounds for total q estimation
#' @param Rel A list of inter-stock relationships see slot Rel of MOM object class
#' @author T.Carruthers
#' @keywords internal
getq_multi_MICE<-function(x,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF, MPA,CatchFrac, bounds= c(1e-05, 15),Rel){

  Nx<-N[x,,,,]
  VFx<-VF[x,,,,]
  FretAx<-FretA[x,,,,]
  NIL(StockPars,"K")

  Kx<-matrix(unlist(lapply(StockPars,function(dat)dat['K'])),ncol=np)[x,]
  Linfx<-matrix(unlist(lapply(StockPars,function(dat)dat['Linf'])),ncol=np)[x,]
  t0x<-matrix(unlist(lapply(StockPars,function(dat)dat['t0'])),ncol=np)[x,]
  Mx<-matrix(unlist(lapply(StockPars,function(dat)dat['M'])),ncol=np)[x,]
  R0x<-matrix(unlist(lapply(StockPars,function(dat)dat['R0'])),ncol=np)[x,]

  hsx<-matrix(unlist(lapply(StockPars,function(dat)dat['hs'])),ncol=np)[x,]
  ax<-matrix(unlist(lapply(StockPars,function(dat)dat['a'])),ncol=np)[1,]
  bx<-matrix(unlist(lapply(StockPars,function(dat)dat['b'])),ncol=np)[1,]
  SRrelx<-matrix(unlist(lapply(StockPars,function(dat)dat['SRrel'])),ncol=np)[x,]

  distx<-SSBpRx<-R0ax<-aRx<-bRx<-array(NA,c(np,nareas))
  Perrx<-array(NA,c(np,nyears+maxage))
  movx<-array(NA,c(np,maxage,nareas,nareas))

  for(p in 1:np){
    distx[p,]<-StockPars[[p]]$R0a[x,]/sum(StockPars[[p]]$R0a[x,])
    Perrx[p,]<-StockPars[[p]]$Perr_y[x,1:(nyears+maxage)]
    movx[p,,,]<-StockPars[[p]]$mov[x,,,]
    SSBpRx[p,]<-StockPars[[p]]$SSBpR[x,]
    R0ax[p,]<-StockPars[[p]]$R0a[x,]
    aRx[p,]<-StockPars[[p]]$aR[x,]
    bRx[p,]<-StockPars[[p]]$bR[x,]
  }

  M_ageArrayx<-Mat_agex<-array(NA,c(np,maxage,nyears))
  Effind<-array(NA,c(np,nf,nyears))
  Spat_targ<-array(NA,c(np,nf))

  for(p in 1:np){
    Mat_agex[p,,]<-StockPars[[p]]$Mat_age[x,,1:nyears]
    M_ageArrayx[p,,]<-StockPars[[p]]$M_ageArray[x,,1:nyears]
    Effind[p,,]<-t(matrix(unlist(lapply(FleetPars[[p]],function(dat,x)dat['Find'][[1]][x,],x=x)),ncol=nf))
    Spat_targ[p,]<-unlist(lapply(FleetPars[[p]],function(dat,x)dat['Spat_targ'][[1]][x],x=x))
  }

  CF<-t(matrix(unlist(lapply(CatchFrac,function(dat)dat[x,])),nrow=nf))
  Fdist<-CF/Effind[,,nyears] # Catch divided by effort (q proxy)
  Fdist<-Fdist/apply(Fdist,1,sum)    # q ratio proxy (real space)
  par<-c(rep(-5,np),log(Fdist[,2:nf]/(1-Fdist[,2:nf]))) # low initial F followed by logit guess at fraction based on Fdist according to catch fraction in recent year

  depc=matrix(unlist(lapply(StockPars,function(dat)dat['D'])),ncol=np)[x,]
  CFc<-array(NA,c(np,nf))
  for(p in 1:np)CFc[p,]=CatchFrac[[p]][x,]

  opt<-optim(par,qestMICE,
             method="L-BFGS-B",lower=c(rep(log(bounds[1]),np),rep(-5,np*(nf-1))),upper=c(rep(log(bounds[2]),np),rep(5,np*(nf-1))),
             depc=depc,CFc=CFc,mode='opt',np=np,nf=nf,nyears=nyears,nareas=nareas,Nx=Nx,VFx=VFx,FretAx=FretAx,
             Effind=Effind,distx=distx,movx=movx,Spat_targ=Spat_targ,M_ageArrayx=M_ageArrayx,Mat_agex=Mat_agex,Kx=Kx,
             Linfx=Linfx,t0x=t0x,Mx=Mx,R0x=R0x,R0ax=R0ax,SSBpRx=SSBpRx,hsx=hsx,ax=ax,bx=bx,aRx=aRx,bRx=bRx,Perrx=Perrx,SRrelx=SRrelx,Rel=Rel,
             control=list(trace=1))

  out<-qestMICE(par=opt$par, depc=depc,CFc=CFc,mode='calc',np=np,nf=nf,nyears=nyears,nareas=nareas,Nx=Nx,VFx=VFx,FretAx=FretAx,
          Effind=Effind,distx=distx,movx=movx,Spat_targ=Spat_targ,M_ageArrayx=M_ageArrayx,Mat_agex=Mat_agex,Kx=Kx,
          Linfx=Linfx,t0x=t0x,Mx=Mx,R0x=R0x,R0ax=R0ax,SSBpRx=SSBpRx,hsx=hsx,aRx=aRx,bRx=bRx,ax=ax,bx=bx,Perrx=Perrx,SRrelx=SRrelx,Rel=Rel)

  return(out)

}


#' Internal function for optimizing catchability (q) for a MICE model
#'
#' Function returns objective function that fits both stock depletion and catch fraction among fleets
#'
#' @param par Integer, the simulation number
#' @param depc Numeric vector, nstock long of specified stock depletion (SSB now / SSB0)
#' @param CFc Matrix [nstock, nfleet], a catch fraction among fleets (sums to 1 for each stock (row))
#' @param mod Character if 'opt' qestMICE returns the objective function otherwise the fitted values in a list
#' @param nf Integer, number of stocks
#' @param nf Integer, number of fleets
#' @param nyears Integer, number of historical years (unfished til today)
#' @param nareas Integer, number of areas (default is 2)
#' @param Nx Array [stock, age, year, area] of stock numbers
#' @param VFx Array [fleet, age, year, area] of the vulnerability curve
#' @param FretAx Array [fleet, age, year, area] of the retention curve
#' @param Effind Array [fleet, year] of effort
#' @param movx Array [stock,age,area,area] of movement transitions
#' @param Spat_targ Matrix [stock, fleet] of spatial targetting parameter (0 evenly spatial distributed, 1 proportional to vulnerable biomass)
#' @param M_ageArrayx Array [stock, age,year] of Natural mortality rate at age
#' @param Mat_agex Array [stock, age, year] of maturity (spawning fraction) age age
#' @param Kx Vector [stock] of von B growth parameter K
#' @param Linf Vector [stock] of von B asymptotic length parameter Linf
#' @param t0 Vector [stock] of von B theoretical age at zero length (t0)
#' @param Mx Vector [stock] mature natural mortality rate
#' @param R0x Vector [stock] unfished recruitment
#' @param R0ax Matrix [stock, area] unfished recruitment by area
#' @param SSBpRx Matrix [stock, area] spawning biomass per recruit by area
#' @param hsx Vector [stock] steepness of the stock recruitment curve
#' @param aRx Vector [stock] stock recruitment parameter alpha (for Ricker curve)
#' @param bRx Vector [stock] stock recruitment parameter beta (for Ricker curve)
#' @param ax Vector [stock] weight-length parameter a W=aL^b
#' @param bx Vector [stock] weight-length parameter b W=aL^b
#' @param Perrx Matrix [stock, year] process error - the lognormal factor for recruitment strength
#' @param SRrelx Integer vector [stock] the form of the stock recruitment relationship (1 = Beverton-Holt, 2= Ricker)
#' @param Rel A list of inter-stock relationships see slot Rel of MOM object class
#' @author T.Carruthers
#' @keywords internal
#' @export
qestMICE<-function(par,depc,CFc,mode='opt',np,nf,nyears,nareas,Nx,VFx,FretAx,Effind,distx,movx,Spat_targ,M_ageArrayx,Mat_agex,
                   Kx,Linfx,t0x,Mx,R0x,R0ax,SSBpRx,hsx,aRx, bRx, ax,bx,Perrx,SRrelx,Rel){

  qs<-exp(par[1:np])
  qlogit<-array(0,c(np,nf))
  qlogit[,2:nf]<-par[(np+1):length(par)]
  qfrac<-exp(qlogit)/apply(exp(qlogit),1,sum)

  HistVars<-popdynMICE(qs,qfrac,np,nf,nyears,maxage,Nx,VFx,FretAx,Effind,movx,Spat_targ,M_ageArrayx,Mat_agex,Kx,Linfx,t0x,Mx,R0x,R0ax,SSBpRx,hsx,aRx, bRx,ax,bx,Perrx,SRrelx,Rel)
  # matplot(t(apply(HistVars$SSBx,c(1,3),sum)))
  # matplot(t(apply(HistVars$Nx,c(1,3),sum)))
  # matplot(t(HistVars$Fy))

  SSBest<-apply(HistVars$SSBx,c(1,3),sum)
  deppred<-SSBest[,nyears]/SSBest[,1]
  Cpred0<-array(NA,c(np,nf,maxage,nareas))
  Cind<-TEG(dim(Cpred0))
  Find<-cbind(Cind[,1:3],nyears,Cind[,4]) # p f age y area
  Bind<-Find[,c(1,3:5)]
  Cpred0[Cind]<-HistVars$Bx[Bind]*(1-exp(-HistVars$FMy[Find]))
  Ctot<-apply(Cpred0,1:2,sum)
  Cpred<-Ctot/apply(Ctot,1,sum)

  depOBJ<-sum((log(depc) - log(deppred))^2)
  cOBJ<-sum(log(CFc[,2:nf]/Cpred[,2:nf])^2)

  if(mode=='opt'){
    return(depOBJ+cOBJ)
  }else{
    return(list(qtot=qs,qfrac=qfrac,CFc=CFc,Cpred=Cpred,depc=depc,deppred=deppred))#,Vulnf=Vulnf,Retf=Retf,MPAf=MPAf))
  }

}



#' optimize for catchability (q)
#'
#' Function optimizes catchability (q, where F=qE) required to get to user-specified stock
#' depletion
#'
#' @param x Integer, the simulation number
#' @param D A numeric vector nsim long of sampled depletion
#' @param SSB0 A numeric vector nsim long of total unfished spawning biomass
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas.
#' Only values from the first year (i.e `N[,,1,]`) are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year
#' @param Mat_age An array (dimensions nsim, maxage, proyears+nyears) with the proportion mature for each age-class
#' @param Asize A matrix (dimensions nsim, nareas) with size of each area
#' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year
#' @param FleetP A list of Fleets containing entries for V, retA, Find and Spat_targ
#' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
#' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
#' @param hs A numeric vector nsim long with the steepness values for each simulation
#' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
#' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
#' @param aR A numeric vector nareas long with the Ricker SRR a values
#' @param bR A numeric vector nareas long with the Ricker SRR b values
#' @param bounds A numeric vector of length 2 with bounds for the optimizer
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param MPA A matrix of spatial closures by year
#' @param useCPP logical - use the CPP code? For testing purposes only
#' @author T.Carruthers and A. Hordyk
#' @keywords internal
#' @export
getq_multi <- function(x, D, SSB0, nareas, maxage, Np, pyears, M_ageArray, Mat_age, Asize, Wt_age,
                  FleetP, Perr, mov, SRrel, hs, R0a, SSBpR, aR, bR,
                  bounds = c(1e-05, 15),maxF, MPAc, CFp, useCPP=TRUE) {

  # Prep multi-fleet arrays
  nf<-length(FleetP)
  ns<-dim(FleetP[[1]]$V)[1]
  ay<-dim(FleetP[[1]]$V)[3]
  Vuln=retA=array(NA,c(nf,maxage,ay))
  ny<-pyears

  Effind<-t(matrix(unlist(lapply(FleetP,function(dat,x)dat['Find'][[1]][x,],x=x)),ncol=nf))
  Spat_targc<-unlist(lapply(FleetP,function(dat,x)dat['Spat_targ'][[1]][x],x=x))

  for(ff in 1:nf){ # kind of lazy but not THAT slow

    Vuln[ff,,]<-FleetP[[ff]]$V[x,,]
    retA[ff,,]<-FleetP[[ff]]$retA[x,,]

  }

  Fdist<-CFp[x,]/Effind[,ny] # Catch divided by effort (q proxy)
  Fdist<-Fdist/sum(Fdist)    # q ratio proxy (real space)
  par<-c(-5,log(Fdist[2:nf]/(1-Fdist[2:nf]))) # low initial F followed by logit guess at fraction based on Fdist according to catch fraction in recent year

  opt <- optim(par,optQ_multi, method="L-BFGS-B",lower=c(log(bounds[1]),rep(-5,nf-1)),upper=c(log(bounds[2]),rep(5,nf-1)),depc=D[x], SSB0c=SSB0[x], nareas=nareas, maxage=maxage, Ncurr=Np[x,,1,],
                  pyears=pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=Vuln, Retc=retA, Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x],
                  Effind=Effind,  Spat_targc=Spat_targc, hc=hs[x], R0c=R0a[x,],
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPAc=MPAc, CFc=CFp[x,],useCPP=useCPP,mode='opt')

  out<-optQ_multi(opt$par,depc=D[x], SSB0c=SSB0[x], nareas=nareas, maxage=maxage, Ncurr=Np[x,,1,],
                  pyears=pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=Vuln, Retc=retA, Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x],
                  Effind=Effind,  Spat_targc=Spat_targc, hc=hs[x], R0c=R0a[x,], SSBpRc=SSBpR[x,],
                  aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPAc=MPAc, CFc=CFp[x,],useCPP=useCPP,mode='out')

  return(out)

}


#' Optimize q for a single simulation
#'
#' @param logQ log q
#' @param depc Depletion value
#' @param SSB0c Unfished spawning biomass
#' @param nareas Number of areas
#' @param maxage Maximum age
#' @param Ncurr Current N-at-age
#' @param pyears Number of years to project population dynamics
#' @param M_age M-at-age
#' @param Asize_c Numeric vector (length nareas) with size of each area
#' @param MatAge Maturity-at-age
#' @param WtAge Weight-at-age
#' @param Vuln Vulnerability-at-age
#' @param Retc Retention-at-age
#' @param Prec Recruitment error by year
#' @param movc movement matrix
#' @param SRrelc SR parameter
#' @param Effind Historical fishing effort
#' @param Spat_targc Spatial targetting
#' @param hc Steepness
#' @param R0c Unfished recruitment by area
#' @param SSBpRc Unfished spawning biomass per recruit by area
#' @param aRc Ricker aR
#' @param bRc Ricker bR
#' @param maxF maximum F
#' @param MPA A matrix of spatial closures by year
#' @param useCPP Logical. Use the CPP code?
#' @param model Character 'opt' is optimization and returns an objective funciton to be minimized, anything else returns fitted values
#' @author A. Hordyk
#' @keywords internal
#' @export
optQ_multi <- function(par, depc, SSB0c, nareas, maxage, Ncurr, pyears, M_age,
                 MatAge, Asize_c, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc,
                 R0c, SSBpRc, aRc, bRc, maxF, MPAc, CFc, useCPP, mode = 'opt') {


    nf<-dim(Effind)[1]
    ny<-dim(Effind)[2]
    ay<-dim(Vuln)[3]

    qtot<-exp(par[1])
    qlogit<-c(0,par[2:nf])
    qfrac<-exp(qlogit)/(sum(exp(qlogit)))
    Effdist<-qfrac*Effind
    Efftot<-apply(Effdist,2,sum)

    MPAind<-TEG(c(nf,ny,nareas))
    MPAtemp<-array(1/nf,dim(MPAc))
    MPAtemp[MPAind]=(MPAc[MPAind]*Effdist[MPAind[,1:2]])/Efftot[MPAind[,2]] # weighted by effort and fleet exposure
    MPAf<-apply(MPAtemp,2:3,sum)

    Vulnf<-Retf<-array(NA,c(maxage,ay))
    Vind<-TEG(c(nf,maxage,ny))
    Vtemp<-Vuln[,,1:ny]
    Vtemp[Vind]<-(Vuln[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
    Vulnf[,1:ny]<-apply(Vtemp,2:3,sum)
    Vulnf[,(ny+1):ay]<-Vulnf[,ny]  # Future vulnerability is the same
    Vulnf<-nlz(Vulnf,2,"max")       # normalize to max 1

    Rtemp<-Retc[,,1:ny]
    Rtemp[Vind]<-(Retc[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
    Retf[,1:ny]<-apply(Rtemp,2:3,sum)
    Retf[,(ny+1):ay]<-Retf[,ny]   # Future retention is the same

    Spat_targf<-sum(apply(Effdist,1,sum)*Spat_targc)/sum(Effdist) # Approximation according to historical F by fleet

    simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                        MatAge, WtAge, Vuln=Vulnf, Retc=Retf, Prec, movc, SRrelc, Effind=Efftot, Spat_targc=Spat_targf, hc,
                        R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=qtot, Fapic=0,
                        maxF=maxF, MPA=MPAf, control=1,  SSB0c=SSB0c)



    # 1:Narray 2:Barray 3:SSNarra 4:SBarray 5:VBarray 6:FMarray 7:FMretarray 8:Zarray
    ssb <- sum(simpop[[4]][,pyears,])
    BB<- apply(simpop[[2]][,pyears,],1,sum)  # age area
    Cf<-array(rep(BB,each=nf)*Vuln[,,ny]*Effdist[,ny],c(nf,maxage))
    Cpred<-apply(Cf,1,sum)/sum(Cf)

    depOBJ<-50*(log(depc) - log(ssb/SSB0c))^2
    cOBJ<-sum(log(Cpred[2:nf]/CFc[2:nf])^2)

    if(mode=='opt'){
      return(depOBJ+cOBJ)
    }else{
      return(list(qtot=qtot,qfrac=qfrac,CFc=CFc,Cpred=Cpred,depc=depc,deppred=ssb/SSB0c,Vulnf=Vulnf,Retf=Retf,MPAf=MPAf))
    }

}



