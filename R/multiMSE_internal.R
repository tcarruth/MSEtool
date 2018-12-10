# multiMSE funcs


#' Slot in list: get the slot values from a list of objects
#'
#' Create of vector of values that correspond with a slot in a list of objects
#'
#' @param listy A list of objects
#' @param sloty A character vector representing the slot name
#' @author T. Carruthers
#' @export
SIL<-function(listy,sloty) {
  if(class(listy[[1]])=="list"){
    out<-NULL
    for(i in 1:length(listy))out<-c(out,unlist(lapply(listy[[i]],function(x)slot(x,sloty))))
  }else{
    out<-unlist(lapply(listy,function(x)slot(x,sloty)))
  }
  out
}


#' Item in list: get the list values from a list of lists
#'
#' Create of vector of values that correspond with a slot in a list of objects
#'
#' @param listy A list of objects
#' @param namey A character vector representing the list item's name
#' @author T. Carruthers
#' @export
NIL<-function(listy,namey){
  unlist(lapply(listy,function(x)x[namey]))
}

#' Sum over list: get the list values from a list of lists
#'
#' Create of vector of values that correspond with a named position in a list of objects
#'
#' @param listy A list of objects
#' @param namey A character vector representing the list item's name
#' @author T. Carruthers
#' @export
SOL<-function(listy,namey){
  out<-listy[[1]][namey][[1]]
  for(i in 2:length(listy))out<-out+listy[[i]][namey][[1]]
  out
}

#' Toms expand grid
#'
#' Create an indexing grid from just a vector of maximum dimension sizes
#'
#' @param vec A vector of maximum array sizes
#' @author T. Carruthers
#' @export
TEG<-function(vec){ # make index for list calculation
  dim<-new('list')
  ndims<-length(vec)
  for(i in 1:ndims)dim[[i]]<-1:vec[i]
  as.matrix(expand.grid(dim))
}



#' Normalize
#'
#' Normalize an array over certain dimensions
#'
#' @param listy A list of objects
#' @param namey A character vector representing the list item's name
#' @author T. Carruthers
#' @export
nlz<-function(arr,dims=NULL,func="max"){
  arrdim<-dim(arr)
  ndim<-length(arrdim)
  if(is.null(dims))dims=1:(ndim-1)
  agg<-apply(arr,dims,func)
  out<-arr
  ind<-TEG(arrdim)
  out[ind]<-out[ind]/agg[ind[,dims]]
  out
}
