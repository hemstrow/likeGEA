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

