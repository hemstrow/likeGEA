#' Simulated example likelihood-GEA data
#'
#' A simulated dataset containing info on 1500 sequenced variants in 200
#' individuals. Contains phred-scaled data, patch info and x/y coords, and
#' environmental variables for each patch.
#'
#'
"example_gl"


#' Example sample metadata
#'
#' A simulated dataset containing patch data, lat/lon coordinates, and
#' environmental data for the 200 individuals in \code{\link{example_gl}}, in
#' an \code{\link[sf]{sf}} format as expected by \code{\link{prepare_metadata}}.
#'
#' Individuals are in rows, sorted identically to genotypic data. Expects
#' columns named:
#' * patch: name of the patch from which an individual sample was collected.
#' * ind_id: individual IDs.
#'
#' All remaining variables must be numeric and are assumed to be uncorrelated
#' environmental variables. The function \code{\link{get_uncorrelated_envs}} can
#' be used to generate uncorrelated variables from a full matrix of envs.
#'
"sample_metadata"
