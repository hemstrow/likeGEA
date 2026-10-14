setup_likeGEA_env <- function(name = "likeGEA", conda_path = "DEFAULT"){

  if(conda_path == "DEFAULT"){
    if(tryCatch(!file.exists(reticulate::conda_binary()), error = function(e) FALSE)){
      reticulate::install_conda()
    }
  } else {
    options(reticulate.conda_binary = conda_path)
  }

  if (reticulate::condaenv_exists(name)) {
    stop("Designated env already exists. Please use a different name or remove.\n")
  } else {
    reticulate::conda_create(name)
    reticulate::conda_install(name, "r-base=4.5.2")

    #========================compilers====================
    reticulate::conda_install(name, packages = c("gcc_linux-64=14.2",
                                                 "gxx_linux-64=14.2",
                                                 "gfortran_linux-64=14.2"))

    # check compilers
    msg <- character()
    check <- reticulate::conda_run2(cmd_line = "R CMD config CC", envname = name, intern = TRUE)
    if(length(grep("conda", check)) == 0){
      msg <- c(msg, "R CMD config CC did not locate conda gnu.")
    }

    check <- reticulate::conda_run2(cmd_line = "R CMD config CXX", envname = name, intern = TRUE)
    if(length(grep("conda", check)) == 0){
      msg <- c(msg, "R CMD config CXX did not locate conda gnu.")
    }

    check <- reticulate::conda_run2(cmd_line = "R CMD config CXX17", envname = name, intern = TRUE)
    if(length(grep("conda", check)) == 0){
      msg <- c(msg, "R CMD config CXX did not locate conda gnu.")
    }

    check <- reticulate::conda_run2(cmd_line = "R CMD config FC", envname = name, intern = TRUE)
    if(length(grep("conda", check)) == 0){
      msg <- c(msg, "R CMD config FC did not locate conda gnu.")
    }

    if(length(msg) > 0){
      stop(paste0(msg, collapse = "\n"))
    }

    reticulate::conda_install(name, channel = "conda-forge", packages = "zlib")

    #==============install cmdstan================
    reticulate::conda_install(name, channel = "conda-forge", packages = "cmdstan")

    #==============install R dependencies==========
    # NOTE: we need to skip the install of rstan for now, since it can cause issues!
    reticulate::conda_install(name, channel = "conda-forge",
                              packages = c("r-cowplot",
                                           "r-geodist",
                                           "r-laplacesdemon",
                                           "r-matrixstats",
                                           "r-posterior",
                                           "r-r.utils",
                                           "r-tictoc",
                                           "r-sf",
                                           "r-terra")) # needed for rnaturalearth later

    # bread, requires direct activation of a R cmd inside the env.
    pkg_inst <- "install.packages(c('bread', 'geodist', 'dplyr', 'gstat', 'stringi', 'rnaturalearth', 'rnaturalearthdata', 'codetools'), repos = 'https://cloud.r-project.org')"
    pkg_inst <- paste("Rscript -e", shQuote(pkg_inst))

    reticulate::conda_run2(
      cmd_line = pkg_inst,
      envname = name)

    # cmdstanr, as before
    cmdstanr_inst <- 'install.packages( "cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))'
    cmdstanr_inst <- paste("Rscript -e", shQuote(cmdstanr_inst))

    reticulate::conda_run2(
      cmd_line = cmdstanr_inst,
      envname = name)


    # check the cmdstan/cmdstanr install
    check <- reticulate::conda_run2(
      cmd_line = paste0("Rscript -e ", shQuote('cmdstanr::cmdstan_version()')),
      envname = name, intern = TRUE)
    check <- check[!grepl("R_HOME", check)]
    if(length(check) > 0){
      cat("Successfully installed cmdstan, version: ", check)
    } else {
      stop("cmdstan version not found by cmdstanr in env.\n")
    }


    check <- reticulate::conda_run2(
      cmd_line = paste0("Rscript -e ", shQuote('cmdstanr::cmdstan_path()')),
      envname = name, intern = TRUE)
    check <- check[!grepl("R_HOME", check)]
    if(length(check) > 0){
      cat("\tpath: ", check, ".\n\t Please ensure that this points to the conda version.\n")
    } else {
      stop("cmdstan not found by cmdstanr in env.\n")
    }

    reticulate::conda_run2(
      cmd_line = paste0("Rscript -e ", shQuote('cmdstanr::check_cmdstan_toolchain()')),
      envname = name)

    #============install likeGEA============
    remotes_inst <- 'install.packages("remotes", repos = c("https://cloud.r-project.org"))'
    remotes_inst <- paste("Rscript -e", shQuote(remotes_inst))

    reticulate::conda_run2(
      cmd_line = remotes_inst,
      envname = name)

    likeGEA_inst <- 'remotes::install_github("hemstrow/likeGEA", dependencies = FALSE)' # need to NOT install rstan.
    likeGEA_inst <- paste("Rscript -e", shQuote(likeGEA_inst))

    reticulate::conda_run2(
      cmd_line = likeGEA_inst,
      envname = name)

    #============install rnaturalearthhires==========
    rnehr_inst <- 'remotes::install_github("ropensci/rnaturalearthhires")' # this one is only on github, installed here since we just did remotes
    rnehr_inst <- paste("Rscript -e", shQuote(rnehr_inst))

    reticulate::conda_run2(
      cmd_line = rnehr_inst,
      envname = name)

    #============install rstan============
    # cmake and the RccParallel via install.packages?
    reticulate::conda_install(name, channel = "conda-forge", packages = "cmake") # needed for RccParallel and StanHeaders


    rstan_inst <- 'install.packages("rstan", repos = c("https://cloud.r-project.org"))'
    rstan_inst <- paste("Rscript -e", shQuote(rstan_inst))

    reticulate::conda_run2(
      cmd_line = rstan_inst,
      envname = name)
  }

  return(name)
}
