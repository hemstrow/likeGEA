#' Format a beagle file as phred-scaled likelihoods.
#'
#' Takes a beagle formatted input file and converts it to a matrix of phred-
#' scaled likelihoods, as needed by \code{likeGEA}. See details.
#'
#' The input file format required by this package has one row for each locus and
#' one column for each individual. Phred-scaled likelihoods (PL) for the three
#' possible genotypes at each locus/individual are comma seperated (major
#' homozygote, heterozygote, minor heterozygote). A higher PL constitutes a
#' lower likelihood; by convention the minimum PL is subtracted from each
#' triplet such that the most likely genotype has a PL of 0.
#'
#' Phred-likelihoods can be calculated according to:
#' \deqn{PL(\text{genotype}) = -10 \times \log_{10}\left( P(\text{Data} \mid \text{Genotype}) \right)}
#'
#' Where \eqn{P(\text{Data} \mid \text{Genotype})} is the a genotype likelihood
#' as produced by something like \code{ANGSD}'s \code{doGlf}.
#'
#' @param file character. File path to the beagle or beagle.gz file containing
#'   the raw data.
#' @param lines numeric, default NULL. Optional vector of numbers noting which
#'   lines in the beagle file to keep. This mut be sequential due to
#'   \code{\link[data.table]{fread}}'s restirctions.
#' @param ind_order numeric, default NULL. Numeric vector indicating final
#'   sorted order to use for individuals the output matrix (\code{c(2,1,3)}, for
#'   example, would put the second individual in \code{file} first, then the
#'   first individual second, and then the third individual last in the final
#'   matrix, for example). If NULL, the order will be unchanged.
#' @param ind_names character, default NULL. Names for individuals in the final
#'   matrix. If \code{ind_order} is provided, this should match the final order
#'   dictated by \code{ind_order}; otherwise, this should match the order
#'   in \code{file}.
#'
#' @return A \eqn{p \times n} matrix, where p is the number of loci and n is
#'   the number of individuals. See details for format.
#'
#' @author William Hemstrom
#'
#' pl <- function("path/to/example/data/once/Nikunj/adds/it.beagle.gz",
#'                lines = 1:2000) # first 2000 lines
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


prepare_pl_batches <- function(file,
                               outdir = "./data",
                               nsnps = 10000,
                               target_order = NULL,
                               individual_IDs = NULL,
                               ind_names = NULL){

  if(!is.null(target_order)){
    if(is.character(target_order)){
      if(!is.null(individual_IDs)){
        stop("Individual_IDs must be provided to re-sort based on IDs.\n")
      }
      if(any(!target_order %in% individual_IDs)){
        stop("All target samples must be present in 'individual_IDs' to re-sort based on IDs.\n")
      }
    }
  }

  realized_order <- readRDS("data/0_ind_IDs.RDS")$ind_id
  selected_snps <- readLines("data/selected_snps.txt")
  selected_snps <- as.numeric(selected_snps)
  if(!all(selected_snps[1]:selected_snps[length(selected_snps)] == selected_snps)){
    stop("selected_snps.txt must give snps as a continuous block for read efficiency.\n")
  }

  resort <- match(target_order, realized_order)
  tpl <- format_pl("data/genotypes.beagle.gz", selected_snps, resort, target_order)
  saveRDS(tpl, "data/1_phred_score.RDS")
  rm(tpl)
}
