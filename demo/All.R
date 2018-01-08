devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("demo", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

source("Medfly.R")

source("KlanjscekDEB.R")

source("MartinsDEB.R")

source("StageStructuredBiomass.R")

source("deRoosPersson.R")

source("Indet_growth.R")

source("deRoosPersson5.R")

source("Indet_growth5.R")

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
