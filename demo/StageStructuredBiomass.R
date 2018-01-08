devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

modelname="StageStructuredBiomass"

cat('\n\n\nDetection of the transcritical bifurcation of the structured population as a function of juvenile mortality\n\n');
cmd = paste0('output1 <- PSPMequi("', modelname, '", "EQ", c(0.8, 30.0), -0.1, c("Muj", 0.1, 1.0), options = c("popZE", "0"), clean = TRUE)')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

plot.new()
par(mar = c(4, 4, 2, 4), tcl=0.4)
plot(1, 1, type="l", xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=c(0.0, 0.6), ylim=c(0.0, 12.5), xlab="", ylab="")
axis(1, at=(0:6)*0.1, labels=TRUE)
axis(2, at=(0:4)*3, labels=TRUE, las=2)
mtext("Juvenile mortality", 1, line=2.5, cex=1.3)
mtext("Juvenile biomass", 2, line=2.5, cex=1.3, col=rgb(0,0,0.6))

lines(output1$curvepoints[,1], output1$curvepoints[,5], type="l", col=rgb(0,0,0.6), lwd=3)
points(output1$bifpoints[,1], output1$bifpoints[,5], col="red", pch=8, lwd=2)
text(1.06*output1$bifpoints[,1], output1$bifpoints[,5], output1$biftype, pos=3, offset=0.45, cex=0.8)

par(new = TRUE)
plot(output1$curvepoints[,1], output1$curvepoints[,6], type="l", col=rgb(0.6,0,0), lwd=3, axes = FALSE, xaxs="i", yaxs="i", bty = "n", xlim=c(0.0, 0.6), ylim=c(0, 4.5), xlab = "", ylab = "")
axis(side=4, at =(0:4), las=2)
mtext("Adult biomass", side=4, line=2.5, cex=1.3, col=rgb(0.6,0,0))
# points(output1$bifpoints[,1], output1$bifpoints[,6], col="red", pch=8, lwd=2)
# text(output1$bifpoints[,1], output1$bifpoints[,6], output1$biftype, pos=1, offset=0.35)

cat("\n\n\nContinuation of the non-trivial equilibrium of the structured population starting from the transcritical bifurcation\n\n")
cmd = paste0('output2 <- PSPMequi("', modelname, '", "EQ", c(output1$bifpoints[c(1, 2, 3)]), -0.5, c("Muj", 0, 2.0))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot.new()
par(mar = c(4, 4, 2, 4), tcl=0.4)
plot(1, 1, type="l", xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=c(0.0, 0.6), ylim=c(0.0, 4.0), xlab="", ylab="")
axis(1, at=(0:7)*0.1, labels=TRUE)
axis(2, at=(0:4), labels=TRUE, las=2)
mtext("Juvenile mortality", 1, line=2.5, cex=1.3)
mtext("Juvenile biomass", 2, line=2.5, cex=1.3, col=rgb(0,0,0.6))

lines(output1$curvepoints[,1], output1$curvepoints[,5], type="l", col=rgb(0,0,0.6), lwd=3)
points(output1$bifpoints[,1], output1$bifpoints[,5], col="red", pch=8, lwd=2)
text(1.06*output1$bifpoints[,1], output1$bifpoints[,5], output1$biftype, pos=3, offset=0.45, cex=0.8)
lines(output2$curvepoints[,1], output2$curvepoints[,5], type="l", col=rgb(0,0,0.6), lwd=3)

par(new = TRUE)
plot(1, 1, type="l", axes = FALSE, xaxs="i", yaxs="i", bty = "n", xlim=c(0.0, 0.6), ylim=c(0, 12.0), xlab = "", ylab = "")
axis(side=4, at =(0:4)*3, las=2)
mtext("Adult biomass", side=4, line=2.5, cex=1.3, col=rgb(0.6,0,0))
lines(output1$curvepoints[,1], output1$curvepoints[,6], type="l", col=rgb(0.6,0,0), lwd=3)
lines(output2$curvepoints[,1], output2$curvepoints[,6], type="l", col=rgb(0.6,0,0), lwd=3)

cat("\n\n\nContinuation of the transcritical bifurcation point of the structured population as a function of juvenile mortality and maxim resource density\n\n")
cmd = paste0('output3 <- PSPMequi("', modelname, '", "BP", c(output1$bifpoints[c(1, 2)], 30.0), -0.5, c("Muj", 0, 2.0, "Rmax", 0, 50), options=c("popBP", "0"))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot.new()
par(mar = c(4, 4, 2, 4), tcl=0.4)
plot(1, 1, type="l", xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=c(0.0, 0.6), ylim=c(0.0, 30.0), xlab="", ylab="")
axis(1, at=(0:7)*0.1, labels=TRUE)
axis(2, at=(0:3)*10, labels=TRUE, las=2)
mtext("Juvenile mortality", 1, line=2.5, cex=1.3)
mtext("Maximum resource density", 2, line=2.5, cex=1.3)
lines(output3$curvepoints[,1], output3$curvepoints[,4], type="l", col=rgb(0.6,0,0), lwd=3)

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
