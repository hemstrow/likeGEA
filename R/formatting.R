#' Convert beagle files to phred-scaled likelihoods.
#' 
#' Takes a "beagle" formatted genotype-likelihood file with three columns per
#' locus (reference homozygote, heterozygote, alternate homozygote) and converts
#' to one-column per locus, phred-scaled, comma-seperated likelihoods.
#' 
#' Phred-scaling converts 
#' 
#' @param file character. Path to a beagle file containing genotype information.
#'   Can be g-zipped.
#' @param lines numeric, default NULL. An optional, sequential vector of rows
#'   to read (ignoring the header). \code{20:30}, for example, would read lines
#'   21 to 31 (loci 20-30 skipping the header). Larger sets of rows slow down
#'   compute times and require substantially more memory.
#' @param ind_order numeric, default NULL. An optional vector giving the final
#'   order of individuals to use in the output. Mostly for internal use.
#' @param ind_names character, default NULL. An optional vector giving the names
#'   of the individuals which will be given as column names. If \code{ind_order}
#'   is also provided, \code{ind_names} should be the final order after any
#'   re-sorting
#'   
#' @returns A data.table containing one row for each locus and one column for
#'   each locus. Genotypes are phred-scaled and comma-seperated. See details.
#'  
#' @author William Hemstrom
#' @export
#' 
#' @examples
#' pl <- function("path/to/example/data/once/Nikunj/adds/it.beagle.gz",
#'                lines = 1:2000, # first 2000 lines
#'                ind_names = c(""))
#' 
beagle2pl <- function(file, lines = NULL, ind_order = NULL, ind_names = NULL){
  
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
  
  if(!is.null(ind_order)){
    dfinal <- dfinal[,ind_order]
  }
  if(!is.null(ind_names)){
    colnames(dfinal) <- ind_names
  }
  
  return(dfinal)
}
