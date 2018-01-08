devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

cmd = 'output1 <- PSPMdemo("KlanjscekDEB", c(0, 1.0, -0.02, 0.4, 1.0), clean=TRUE, force=TRUE)'
str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))
plot(output1$curvepoints[,1], output1$curvepoints[,2], type='l', col=rgb(.6,0,0), lwd=3, xlim=c(0.4,1.0), ylim=c(0.0,0.7), xlab="Food density", ylab="Population growth rate", font.lab=2)

cmd = 'output2 <- PSPMdemo("KlanjscekDEB2", c(0, 1.0, -0.02, 0.4, 1.0), clean=TRUE, force=TRUE)'
str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))
lines(output2$curvepoints[,1], output2$curvepoints[,2], type='l', col=rgb(0,.6,0), lwd=3)

cmd = 'output3 <- PSPMdemo("KlanjscekDEBpulsed", c(0, 1.0, -0.02, 0.4, 1.0), clean=TRUE, force=TRUE)'
str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))
lines(output3$curvepoints[,1], output3$curvepoints[,2], type='l', col=rgb(0,0,.6), lwd=3)

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
