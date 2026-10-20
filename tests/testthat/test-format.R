test_that("get_chunks",{

  # correct on linux with gz
  beagle <- system.file("extdata", "example.beagle.gz", package = "likeGEA")
  chunk_files <- get_chunks(beagle, 100)

  ## ref
  ref <- fread(beagle)
  ## chunks
  chunks <- lapply(chunk_files, fread)
  ## colnames
  for(x in 1:length(chunks)){
    expect_identical(colnames(chunks[[x]]), colnames(ref))
  }
  ## contents
  interleaved <- rbindlist(chunks, idcol = "chunk")[, row_num := seq_len(.N), by = chunk]
  interleaved <- interleaved[order(row_num, chunk)]
  interleaved$chunk <- interleaved$row_num <- NULL
  expect_identical(interleaved, ref)
  file.remove(chunk_files)

  # correct on linux with non-gz
  beagle <- system.file("extdata", "example.beagle.gz", package = "likeGEA")
  file.copy(beagle, "test.beagle.gz")
  system("gunzip test.beagle.gz")
  beagle <- "test.beagle"
  chunk_files <- get_chunks(beagle, 100)

  ## ref
  ref <- fread(beagle)
  ## chunks
  chunks <- lapply(chunk_files, fread)
  ## colnames
  for(x in 1:length(chunks)){
    expect_identical(colnames(chunks[[x]]), colnames(ref))
  }
  ## contents
  interleaved <- rbindlist(chunks, idcol = "chunk")[, row_num := seq_len(.N), by = chunk]
  interleaved <- interleaved[order(row_num, chunk)]
  interleaved$chunk <- interleaved$row_num <- NULL
  expect_identical(interleaved, ref)
  file.remove(chunk_files)
  file.remove("test.beagle")
  if(file.exists("test.beagle.gz")){file.remove("test.beagle.gz")}
})
