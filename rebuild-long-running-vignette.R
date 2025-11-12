old_wd <- getwd()

setwd("vignettes/")
knitr::knit("mycoplasma.Rmd.orig", output = "mycoplasma.Rmd")
knitr::purl("mycoplasma.Rmd.orig", output = "mycoplasma.R")

setwd(old_wd)
