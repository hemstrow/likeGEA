# checks if cmdstanr can find cmdstan and stops if not. Otherwise returns path.
.check_cmdstan <- function(){
  path <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  if(methods::is(path, "try-error")){

    say <- "No cmdstan install found.\n"

    if(!interactive()){
      stop(say)
    }

    cat(say, "")
    cat("Install? (y/n)\n")

    resp <- readLines(n = 1)
    resp <- tolower(resp)
    if(resp == "yes") resp <- "y"
    if(resp != "y"){
      stop(say)
    } else{
      cmdstanr::install_cmdstan()
    }

    path <- try(cmdstanr::cmdstan_version(), silent = TRUE)
    if(methods::is(path, "try-error")){
      stop("cmdstan still cannot be located. Likely installation failure. Try a manual install?")
    }
  }
  return(path)
}

.console_hline <- function(char = "="){
  return(paste0(rep(char, getOption("width")), collapse = ""))
}


.resort_phred <- function(file,
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
