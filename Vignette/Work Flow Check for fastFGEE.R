library(devtools)
#path <- "/Users/loewingergc/Desktop/NIMH Research/git_repos/fastFMM-main"
path <- "/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE"
# remove old study strap package
detach("package:fastFGEE", unload=TRUE)
remove.packages("fastFGEE")

setwd(path)

# document and build package
devtools::document()
devtools::build()

# install new one
install.packages("/Users/loewingergc/Desktop/NIMH Research/git_repos/fastFMM_0.3.0.tar.gz", 
                    repos = NULL, type = "source")
library(fastFMM)

# document and build package -- Need to repeat this again because need smtl_setup run before vignette can be 
# properly built
setwd(path)

devtools::document()
devtools::build()

devtools::check()

# STOP HERE !!!!

