beagle2pl <- function(file, lines = NULL, ind.order = NULL, ind.names = NULL){
  
  if(!is.null(lines)){
    d <- fread(file, skip = lines[1] - 1, nrows = length(lines))
  }
  else{
    d <- fread(file)
  }
  
  header <- d[,1:3]
  d <- d[,-c(1:3)]
  d <- as.matrix(d)
  d[d == 0] <- 0.00000000001
  d <- -10*log10(d)
  d <- round(d)
  storage.mode(d) <- "integer"
  
  # loop through inds to standardize, should be quickish
  inds <- 1:ncol(d)
  inds <- split(inds, rep(1:(ncol(d)/3), each = 3))
  for(i in 1:length(inds)){
    m <- matrixStats::rowMins(d[,inds[[i]]])
    d[,inds[[i]]] <- d[,inds[[i]]] - m
  }
  
  # convert to proper strings
  storage.mode(d) <- "character"
  dfinal <- matrix(as.character(NA), nrow(d), ncol(d)/3)
  for(i in 1:length(inds)){
    dfinal[,i] <- paste0(d[,inds[[i]][1]], ",", d[,inds[[i]][2]], ",", d[,inds[[i]][3]])
  }
  
  dfinal[dfinal == "0,0,0"] <- "."
  
  if(!is.null(ind.order)){
    dfinal <- dfinal[,ind.order]
  }
  if(!is.null(ind.names)){
    colnames(dfinal) <- ind.names
  }
  
  return(dfinal)
}
