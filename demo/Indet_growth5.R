devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

modelname="Indet_growth_5bs.h"

cat('\n\n\nDetection of the ESS value of the bifurcation parameter of the structured population (ingestion exponent)\n\n');
cmd = paste0('output1 <- PSPMequi("', modelname, '", "EQ", c(1.0, 0.22, 0.0), -0.1, c(6, 0.5, 2.0), options=c("popEVO", "0"), clean = TRUE)')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

plot.new()
par(mar = c(4, 5, 2, 5), tcl=0.4)
plot(1, 1, type="l", xaxt="n", yaxt="n", xaxs="i", yaxs="i", log="y", xlim=c(0.5, 1.0), ylim=c(0.2, 0.25), xlab="", ylab="")
axis(1, at=(5:10)*0.1, labels=TRUE)
axis(2, at=0.2+(0:5)*0.01, labels=TRUE, las=2)
mtext("Ingestion exponent", 1, line=2.5, cex=1.3)
mtext("Equilibrium resource density", 2, line=3.5, cex=1.3)

lines(output1$curvepoints[,1], output1$curvepoints[,2], type="l", col=rgb(0,.6,0), lwd=3)
points(output1$bifpoints[,1], output1$bifpoints[,2], col="red", pch=8, lwd=2)
text(output1$bifpoints[,1], output1$bifpoints[,2], output1$biftype, pos=1, offset=0.35)

par(new = TRUE)
plot(output1$curvepoints[,1], output1$curvepoints[,5]+output1$curvepoints[,6], type="l", col="blue", lwd=3, axes = FALSE, xaxs="i", yaxs="i", bty = "n", xlim=c(0.5, 1.0), ylim=c(0.78,0.82), xlab = "", ylab = "")
axis(side=4, at = 0.78+(0:4)*0.01, las=2)
mtext("Consumer biomass", side=4, line=3.5, cex=1.3)
points(output1$bifpoints[,1], output1$bifpoints[,5]+output1$bifpoints[,6], col="red", pch=8, lwd=2)
text(output1$bifpoints[,1], output1$bifpoints[,5]+output1$bifpoints[,6], output1$biftype, pos=1, offset=0.35)

cat("\n\n\nContinuation of the ESS value of the ingestion exponentof the structured population as a function of the first bifurcation parameter (maintenance exponent)\n\n")
cmd = c(paste0('output2a <- PSPMequi("', modelname, '", "ESS", c(1.0, output1$bifpoints[c(2, 3, 1)]), -0.1, c(9, 0.5, 2.0, 6, 0.5, 2.0), options=c("popEVO", "0"))'), paste0('output2b <- PSPMequi("', modelname, '", "ESS", c(1.0, output1$bifpoints[c(2, 3, 1)]), 0.1, c(9, 0.5, 2.0, 6, 0.5, 2.0), options=c("popEVO", "0"))'))

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot.new()
plot(output2a$curvepoints[,1], output2a$curvepoints[,4], type='l', col=rgb(0.6,0,0), lwd=3, xlim=c(0.5,2.0), ylim=c(0.5,2.0), xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="", font.lab=2)
lines(output2b$curvepoints[,1], output2b$curvepoints[,4], type='l', col=rgb(0.6,0,0), lwd=3)
axis(1, at=0.5+(0:4)*0.5, labels=TRUE)
axis(2, at=0.5+(0:4)*0.5, labels=TRUE, las=2)
mtext("Maintenance exponent", 1, line=2.5, cex=1.3)
mtext("Ingestion exponent", 2, line=3.5, cex=1.3)

cat("\n\n\nConstruction of the PIP of the resident and mutant value of the first parameter of the structured population (ingestion exponent)\n\n")
cmd = c(paste0('output3 <- PSPMequi("', modelname, '", "PIP", c(output1$bifpoints[c(1, 2, 3, 1)]), 0.1, c(6, 0.5, 2.0, 6, 0.5, 2.0), options=c("popEVO", "0"))'), paste0('output4 <- PSPMequi("', modelname, '", "PIP", c(output1$bifpoints[c(1, 2, 3, 1)]), -0.1, c(6, 0.5, 2.0, 6, 0.5, 2.0), options=c("popEVO", "0"))'))


str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot(output3$curvepoints[,1], output3$curvepoints[,4], type='l', col=rgb(0.6,0,0), lwd=3, xlim=c(0.5,1.5), ylim=c(0.5,1.5), xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="", font.lab=2)
axis(1, at=0.6+(0:4)*0.2, labels=TRUE)
axis(2, at=0.6+(0:4)*0.2, labels=TRUE, las=2)
mtext("Resident ingestion exponent", 1, line=2.5, cex=1.3)
mtext("Mutant ingestion exponent", 2, line=3.5, cex=1.3)
lines(output4$curvepoints[,1], output4$curvepoints[,4], type='l', col=rgb(0.6,0,0), lwd=3)
lines(c(0.5,1.5), c(0.5,1.5), type='l', col="black", lwd=1)

cat("\n\n\nSimulating the evolution in the ingestion exponent over evolutionary time\n\n")
cmd = paste0('output6 <- PSPMevodyn("', modelname, '", c(0.22, 0.03554, 1.0), c(0.05, 10), c(6, 0.5, 1.5), options=c("popEVO", "0"))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot.new()
plot(output6$curvepoints[,1], output6$curvepoints[,4], type='l', col=rgb(0.6,0,0), lwd=3, xlab="", ylab="", font.lab=2)
mtext("Evolutionary time", 1, line=2.5, cex=1.3)
mtext("Ingestion exponent", 2, line=3.5, cex=1.3)

cat("\n\n\nSimulating the simultaneous evolution of the ingestion and maintenance exponent over evolutionary time\n\n")
cmd = paste0('output7 <- PSPMevodyn("', modelname, '", c(0.22, 0.03554, 1.0, 1.0), c(0.05, 100), c(6, 0.5, 1.5, 9, 0.5, 1.5), options=c("popEVO", "0"))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot.new()
plot(output7$curvepoints[,1], output7$curvepoints[,4], type='l', col=rgb(0.6,0,0), lwd=3, xlab="", ylab="", font.lab=2, xlim=c(0.0,60.0), ylim=c(0.5,1.1))
lines(output7$curvepoints[,1], output7$curvepoints[,5], type='l', col=rgb(0,0.6,0), lwd=3)
mtext("Evolutionary time", 1, line=2.5, cex=1.3)
mtext("Ingestion/Maintenance exponent", 2, line=3.5, cex=1.3)

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
