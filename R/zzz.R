.onAttach <- function(libname, pkgname){
  version <- utils::packageVersion("likeGEA")
  packageStartupMessage("likeGEA version: ", version, " loaded.\n", .console_hline(), "\n")
  if(grepl("\\.9[0-9]{3}$", version)){
    packageStartupMessage("This is a development build of likeGEA.\nWhile this version is passing all tests, there may still be bugs present.\nPlease report any issues to the github issues page (https://github.com/hemstrow/likeGEA/issues)!\n\n\n")
    packageStartupMessage("Check out the NEWS page on github for update notes (https://github.com/hemstrow/likeGEA/blob/dev/NEWS.md)\n", .console_hline(), "\n")
  }
  else{
    packageStartupMessage("Please report any issues to the github issues page (https://github.com/hemstrow/likeGEA/issues)!\n\n\n")
    packageStartupMessage("Check out the NEWS page on github for update notes (https://github.com/hemstrow/likeGEA/blob/master/NEWS.md)\n", .console_hline(), "\n")

  }
}
