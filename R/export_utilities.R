#' Generate an example SLURM array job script and handlers to run
#' \code{\link{likeGEA}} on a SLURM HPC.
#'
#' Please note: the resulting shell script, located in the \code{output_dir}
#' named \code{example_slurm_array_job.sh} should be edited to fit the user's
#' SLURM configuration. See notes at the top of the file!
#'
#' @param beagle_files vector of file paths to \code{.beagle}/\code{.beagle.gz}
#'   files to run. These will be formatted to phred files for running by the
#'   generated SLURM script for more effective parallelization.
#' @param sample_metadata either a named list of data.framescontaining the
#'   needed metadata in the correct format expected by \code{\link{likeGEA}}
#'   (see \code{\link{example_gl}}, although the \code{$likelihoods} object is
#'   not needed.) or an input \code{\link[sf]{sf}} object like
#'   \code{\link{sample_metadata}}.
#' @param output_dir character, default ".". A path to the director where final
#'   \code{\link{likeGEA}} outputs will be stored in named sub-directories for
#'   each file in \code{beagle_files}.
#' @param conda_env character, default likeGEA. A named \code{conda} environment
#'   to load at the start of the script. Note that if a conda module etc. needs
#'   to be loaded for this env to load, the resulting SLURM array-job shell
#'   script must be modified. See \code{\link{setup_likeGEA_env}} for help
#'   setting up a working \code{conda} environment.
#' @param control_file_prefix character, default "./control/". A prefix indicating
#'   where the control files (including meta-data and the shell script) will be
#'   located. Either a directory ending in "/" or a file basename prefix.
#' @param env_mean_method character, default "unweighted". Weighting scheme for
#'   mean envs if a \code{\link{sf}} object is provided to \code{sample_metadata}.
#'   See \code{\link{prepare_metadata}} for details and options.
#' @param ... additional arguments passed to \code{\link{likeGEA}}.
#'
#' @export
#' @author William Hemstrom
example_slurm_job <- function(beagle_files,
                              sample_metadata,
                              output_dir = ".",
                              conda_env = "likeGEA",
                              control_file_prefix = "./control/",
                              env_mean_method = "unweighted",
                              ...){

  browser()
  beagle_files <- normalizePath(beagle_files)
  bad_beagles <- !file.exists(beagle_files)
  if(any(bad_beagles)){
    stop(paste0("Some beagle files not found:\n", paste0(bad_beagles, collapse = "\n\n")))
  }

  nfiles <- length(beagle_files)

  # create and sanitize output dir
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  # sanitize control file prefix (expects trailing / IF A DIRECTORY)
  add_trail <-  grepl("[/\\\\]$", control_file_prefix)
  control_file_prefix <- normalizePath(control_file_prefix, mustWork = FALSE)
  if(add_trail) control_file_prefix <- paste0(control_file_prefix, .Platform$file.sep)

  ## make base dir if needed
  bd <- dirname(control_file_prefix)
  if(!dir.exists(bd)) dir.create(bd)

  # if a dir, make it if needed
  if(grepl("[/\\\\]$", control_file_prefix)){
    if(!dir.exists(control_file_prefix)){
      dir.create(control_file_prefix)
    }
  }

  # save a beagle list
  bfp <- normalizePath(paste0(control_file_prefix, "beaglelist.txt"), mustWork = FALSE)
  writeLines(beagle_files, bfp, sep = "\n")
  beagle_files <- bfp

  # format and save the metadata
  if(is.list(sample_metadata) & !methods::is(sample_metadata, "sf")){
    if(all(c("patch_info", "patches_xy", "env") %in% names(sample_metadata))){

      smpath <- normalizePath(paste0(control_file_prefix, "sample_metadata.RDS"), mustWork = FALSE)
      saveRDS(sample_metadata[c("patch_info", "patches_xy", "env")], smpath)
      sample_metadata <- smpath

    } else{
      stop(paste0("'sample_metadata' is a list, but does not contain needed parts:",
                  paste0(c("patch_info", "patches_xy", "env"), collapse = ", "), "\n"))
    }
  } else {
    sample_metadata <- prepare_metadata(sample_metadata, env_mean_method)
    smpath <- normalizePath(paste0(control_file_prefix, "sample_metadata.RDS"), mustWork = FALSE)
    saveRDS(sample_metadata, smpath)
    sample_metadata <- smpath
  }

  # read in the shell script
  browser()
  path <- system.file("extdata", "example_slurm_array_job.sh", package = "likeGEA")
  sh <- paste0(control_file_prefix, "example_slurm_array_job.sh")
  sh <- normalizePath(sh, mustWork = FALSE)
  file.copy(path, sh, overwrite = TRUE)


  # generate the parameter file
  args <- list(...)
  control_parmfile <- paste0(control_file_prefix, "parmfile.RDS")
  normalizePath(control_parmfile, mustWork = FALSE)
  saveRDS(args, control_parmfile)

  # copy over the Rscript
  path <- system.file("extdata", "example_array_rscript.R", package = "likeGEA")
  control_rscript <- paste0(control_file_prefix, "example_array_rscript.R")
  control_rscript <- normalizePath(control_rscript, mustWork = FALSE)
  file.copy(path, control_rscript, overwrite = TRUE)


  # substitute in the provided data
  browser()
  xfun::gsub_file(sh, "BEAG_REPLACE", beagle_files)
  xfun::gsub_file(sh, "META_REPLACE", sample_metadata)
  xfun::gsub_file(sh, "PARAM_REPLACE", control_parmfile)
  xfun::gsub_file(sh, "DIR_REPLACE", output_dir)
  xfun::gsub_file(sh, "ENV_REPLACE", conda_env)
  xfun::gsub_file(sh, "RSCRIPT_REPLACE", control_rscript)
  xfun::gsub_file(sh, "ARRAY_REPLACE", as.character(nfiles))

  return(sh)
}
