likelihoods <- readRDS("/mnt/c/Users/hemst/Downloads/R_package/R_package/example/gl/data/1_phred_score.RDS")
patches <- readRDS("/mnt/c/Users/hemst/Downloads/R_package/R_package/example/gl/data/2_ind_patches.RDS")
patches_xy <- readRDS("/mnt/c/Users/hemst/Downloads/R_package/R_package/example/common/2_patches_xy.RDS")
env <- readRDS("/mnt/c/Users/hemst/Downloads/R_package/R_package/example/common/3_env.RDS")

example_gl <- list(likelihoods = likelihoods,
                   patch_info = patches,
                   patches_xy = patches_xy,
                   env = env)
usethis::use_data(example_gl)
