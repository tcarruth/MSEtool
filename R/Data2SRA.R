

# Needed: a function that converts a DLMtool data object to SRAscope inputs
Data2SRA<-function(dat,OM){ # dat is DLMtool class Data
  
  data<-Process_nonindices() # not yet completed
  data<-c(data,Process_Indices(dat,OM) # Code to get other SRA data from dat object into a list object for SRAscope
  data # 
  
  # for use in (for example): 
  # out<-SRA_scope(OM = OM,
    #             data = data,
    #             mean_fit = TRUE,
    #             cores = ncpus,
    #             plusgroup = T,
    #             OMeff = OMeff,
    #             ESS = data$ESS,
    #             s_selectivity = data$s_selectivity, 
    #             selectivity = data$selectivity,
    #             s_vul_par = data$s_vul_par, s_map_vul_par = data$s_map_vul_par,
    #             vul_par = data$vul_par, map_vul_par = data$map_vul_par,
    #             LWT = data$LWT,
    #             max_F=data$max_F,
    #             control=list(eval.max=1E4, iter.max=1E4, abs.tol=1e-6))
  
}

# experimental code to format index data for a stochastic SRA run where the originate from a DLMtool:Data object 

Process_Indices<-function(dat,OM){  
  
  # https://cran.r-project.org/web/packages/MSEtool/vignettes/SRA_scope_sel.html
  # Index stacking ---------------------------
  Index<-NULL
  I_sd<-NULL
  Itype<-NULL
  sel_block<-NULL
  vul_par<-NULL
  map_vul_par<-NULL
  s_vul_par<-NULL
  s_map_vul_par<-NULL
  s_selectivity<-NULL
  ny<-length(dat@Year)
  fleetno<-0
  
  # Total biomass index -------------
  temp<-slot(dat,"Ind")
  if(!all(is.na(temp))){
    Index<-rbind(Index,temp)
    temp<-slot(dat,"CV_Ind")
    if(!all(is.na(temp))){
      temp[is.na(temp)]<-mean(temp,na.rm=T)
      I_sd<-rbind(I_sd,temp)
    }else{
      I_sd<-rbind(I_sd,rep(mean(OM@Iobs),length(temp)))
    }    
    Itype<-c(Itype,"B")
    
    #selectivity
    temp<-rep(NA,dat@MaxAge)
    s_vul_par<-cbind(s_vul_par,temp)
    
    sel_block<-cbind(sel_block,rep(fleetno,ny))
    
    s_selectivity<-c(s_selectivity,"dome")
    
  }
  
  # Total spawning biomass index ---
  temp<-slot(dat,"SpInd")
  if(!all(is.na(temp))){
    Index<-rbind(Index,temp)
    temp<-slot(dat,"CV_SpInd")
    if(!all(is.na(temp))){
      I_sd<-rbind(I_sd,temp)
    }else{
      I_sd<-rbind(I_sd,rep(mean(OM@Iobs),length(temp)))
    }    
    Itype<-c(Itype,"SSB")
    
    # selectivity
    temp<-rep(NA,dat@MaxAge)
    s_vul_par<-cbind(s_vul_par,temp)
    
    sel_block<-cbind(sel_block,rep(fleetno,ny))
    s_selectivity<-c(s_selectivity,"log")
  }
  
  # Vulnerable biomass (according to pars in the data sheet)
  
  temp<-slot(dat,"VInd")
  if(!all(is.na(temp))){
    
    ny<-ncol(temp)
    Index<-rbind(Index,temp)
    temp<-slot(dat,"CV_VInd")
    
    if(!all(is.na(temp))){
      I_sd<-rbind(I_sd,temp)
    }else{
      I_sd<-rbind(I_sd,rep(mean(OM@Iobs),length(temp)))
    }    
    fleetno<-fleetno+1
    Itype<-c(Itype,fleetno)
    
    # selectivity
    temp<-c(mean(OM@LFS),mean(OM@L5),mean(OM@Vmaxlen))
    vul_par<-cbind(vul_par,temp)
    
    temp<-rep(NA,dat@MaxAge)
    s_vul_par<-cbind(s_vul_par,temp)
    sel_block<-cbind(sel_block,rep(fleetno,ny))
    s_selectivity<-c(s_selectivity,"dome")
    
  }
  
  temp<-slot(dat,"AddInd")[1,,]
  if(!all(is.na(temp))){
    ny<-ncol(temp)
    nadd<-nrow(temp)
    Index<-rbind(Index,temp)
    temp<-slot(dat,"CV_AddInd")[1,,]
    I_sd<-rbind(I_sd,temp)
    temp<-slot(dat,"AddIndV")[1,,]
    s_vul_par<-cbind(s_vul_par,t(temp))
    temp<-slot(dat,"AddIndType")
    
    if(is.na(temp))temp<-rep(3,nadd) # default to 3 (vulnerable index type)
    for(i in 1:length(temp)){
      #fleetno<-fleetno+1
      s_selectivity<-c(s_selectivity,"free")
      #sel_block<-cbind(sel_block,rep(fleetno,ny))
      if(temp[i]==1){
        Itype<-c(Itype,"B")
      }else if(temp[i]==2){
        Itype<-c(Itype,"SSB")
      }else{
        Itype<-c(Itype,'est')
      }
    } 
  }
  
  if(all(is.na(s_vul_par)))s_vul_par=NULL
  if(all(is.na(vul_par)))vul_par=NULL
  
  if(!is.null(s_vul_par))s_map_vul_par=array(NA,dim(s_vul_par))
  if(!is.null(vul_par))map_vul_par=array(NA,dim(vul_par))
  
  #if(all(PanelState[[1]][[10]]==c(T,F,F,F))){
  #  selectivity='logistic'
  #  AM("Conditioning operating model estimating logistic ('flat-topped') selectivity based on Fishery Question 11")
  #}else{  
  #  selectivity='dome'
  #  AM("Conditioning operating model allowing for the estimation of dome shaped selectivity based on Fishery Question 11")
  #}  
  
  list(Index=t(Index), I_sd=t(I_sd), I_type=Itype, 
       s_vul_par=s_vul_par, s_map_vul_par=s_map_vul_par,
       vul_par=vul_par, map_vul_par=map_vul_par,
       #nsel_block=ncol(sel_block), sel_block=sel_block, 
       s_selectivity = s_selectivity, selectivity=selectivity)
  
}
