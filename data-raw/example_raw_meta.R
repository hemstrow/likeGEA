patches <- example_gl$patch_info
patches_xy <- example_gl$patches_xy
env <- example_gl$env
patches_as_provided <- merge(patches, patches_xy, by = "patch", sort = FALSE)
patches_as_provided <- merge(patches_as_provided, env, by = "patch", sort = FALSE)
sample_metadata <- sf::st_as_sf(patches_as_provided, coords = c("lon", "lat"))
sf::st_crs(sample_metadata) <- "EPSG:4326"
usethis::use_data(sample_metadata, overwrite = TRUE)
