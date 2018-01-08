devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

cmd = 'output <- PSPMdemo("Medfly_periodic", c(7, 11.6, 0.1, 11, 20))'
str = readline(paste0("\n> ", cmd, " [y(es)/s(kip)] : "))
if ((str == '') || (str == 'y') || (str == 'Y')) {
  eval(parse(text=cmd))
  plot(output$curvepoints[,1], output$curvepoints[,2], type='l', lwd=2, xlab="Forcing period", ylab="Population growth rate", font.lab=2)
}

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
