
  
SRA_scope_data <- function (OM, Data, ...) {
 
  # < code for making Chist, Index, ML, CAA, CAL, I_sd matrices out of list of Data objects >
  nyears<-OM@nyears
  CYcond<-length(Data@Cat[1,]) != nyears
  if(CYcond) message(paste0("Catch data for a different duration than operating model. Data@Cat[1,] is of length ",length(Data@Cat[1,]),", but OM@nyears is ",nyears))
  if(CYcond)stop("OM, Data, incompatibility")
 
  Chist<-Data@Cat[1, ]
  Index<-Data@Ind[1, ]
  I_sd<-Data@CV_Ind[1, ]
  CAA<-Data@CAA[1, , ]
  CAA<-NULL
  CAL<-Data@CAL[1, , ]
  CAL<-NULL
  ML<-Data@ML[1, ]
  ML<-NULL
  length_bin<-Data@CAL_mids
 
  out<-SRA_scope(OM, Chist=Chist, Index = Index, I_sd = I_sd, CAA = CAA, CAL = CAL, 
                            ML = ML, length_bin = length_bin, ...)
  return(out) 
}


