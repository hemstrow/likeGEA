gl <- example_gl$likelihoods

fl <- phred_to_loglik(gl)

fl1 <- fl[,,1]
fl2 <- fl[,,2]
fl3 <- fl[,,3]


fl <- cbind(fl1, fl2, fl3)
k <- ncol(fl1)
fl <- fl[,order(rep(1:k, 3))]
fl <- as.data.table(fl)
cn <- rep(1:ncol(gl), each = 3)
cn <- paste0("Ind", formatC(cn, width = 3, format = "d", flag = "0"))
colnames(fl) <- cn
snpinfo <- data.table(marker = paste0("CHR1_", 1:nrow(fl)),
                      allele1 = 1,
                      allele2 = 2)

beagle <- cbind(snpinfo, fl)
fwrite(beagle, "inst/extdata/example.beagle.gz", sep = "\t", eol = "\n", row.names = FALSE, col.names = TRUE)
