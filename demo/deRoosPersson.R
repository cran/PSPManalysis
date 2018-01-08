devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

# The function to plot the bifurcation plots in 3 stages
bifplot <- function(stage = 1) {
  devAskNewPage(ask = FALSE)
  graphics.off()
  plot.new()
  par(mfrow = c(3, 1))
  par(cex = 1.0, cex.lab=1.3)
  par(oma = c(0, 0, 0.5, 1.0))
  par(tcl = 0.4)
  par(mgp = c(2, 0.6, 0))

  par(mfg=c(3,1), mar = c(3.5, 4, 0, 0.5))
  plot(1, 1, type="l", xaxt="n", xlim=c(0,0.0004), xaxs="i",yaxs="i", yaxt="n", log="y", ylim=c(5.0E-7, 5.0E-4), xlab="Maximum resource density", ylab="")
  axis(1, at=(0:4)*0.0001, labels=c("0.0", "0.0001", "0.0002", "0.0003", "0.0004"))
  axis(2, at=10^(-6:-4), labels=c(expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2)
  mtext("Resource biomass", 2, line=2.5, cex=1.3)

  lines(output1$curvepoints[,1], output1$curvepoints[,2], type="l", col=rgb(0,.6,0), lwd=3)
  if (stage > 1) {
    lines(output2$curvepoints[,1], output2$curvepoints[,2], type="l", col=rgb(0,.6,0), lwd=3)
    if (stage > 2)
      lines(output3$curvepoints[,1], output3$curvepoints[,2], type="l", col=rgb(0,.6,0), lwd=3)
  }

  points(output1$bifpoints[,1], output1$bifpoints[,2], col="red", pch=8, lwd=2)
  text(output1$bifpoints[,1], output1$bifpoints[,2]-3.5E-6, output1$biftype, pos=4, offset=0.2)
  if (stage > 1) {
    points(output2$bifpoints[,1], output2$bifpoints[,2], col="red", pch=8, lwd=2)
    text(output2$bifpoints[,1], output2$bifpoints[,2], output2$biftype, pos=1, offset=0.5)
    if (stage > 2) {
      points(output3$bifpoints[,1], output3$bifpoints[,2], col="red", pch=8, lwd=2)
      text(output3$bifpoints[,1], output3$bifpoints[,2], output3$biftype, pos=2, offset=0.3)
    }
  }

  par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
  plot(1, 1, type="l", col=rgb(0,0,.6), lwd=3, log="y", xaxt="n", xlim=c(0,0.0004), xlab="", xaxs="i",yaxs="i", yaxt="n", ylim=c(5.0E-7, 1.0E-3), ylab="")
  axis(1, at=(0:4)*0.0001, labels=c("", "", "", "", ""))
  axis(2, at=10^(-7:-4), labels=c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2)
  mtext("Consumer biomass", 2, line=2.5, cex=1.3)

  if (stage > 1) {
    lines(output2$curvepoints[,1], output2$curvepoints[,7]+output2$curvepoints[,8], type="l", col=rgb(0,0,.6), lwd=3)
    lines(output2$curvepoints[,1], output2$curvepoints[,9], type="l", col=rgb(.6,0,0), lwd=3)
    if (stage > 2) {
      lines(output3$curvepoints[,1], output3$curvepoints[,7]+output3$curvepoints[,8], type="l", col=rgb(0,0,.6), lwd=3)
      lines(output3$curvepoints[,1], output3$curvepoints[,9], type="l", col=rgb(.6,0,0), lwd=3)
    }
    points(output2$bifpoints[,1], output2$bifpoints[,7]+output2$bifpoints[,8], col="red", pch=8, lwd=2)
    text(output2$bifpoints[,1], output2$bifpoints[,7]+output2$bifpoints[,8], output2$biftype, pos=3, offset=0.5)

    points(output2$bifpoints[,1], output2$bifpoints[,9], col="red", pch=8, lwd=2)
    text(output2$bifpoints[,1], output2$bifpoints[,9], output2$biftype, pos=3, offset=0.5)
    if (stage > 2) {
      points(output3$bifpoints[,1], output3$bifpoints[,7]+output3$bifpoints[,8], col="red", pch=8, lwd=2)
      text(output3$bifpoints[,1], output3$bifpoints[,7]+output3$bifpoints[,8], output3$biftype, pos=2, offset=0.3)

      points(output3$bifpoints[,1], output3$bifpoints[,9], col="red", pch=8, lwd=2)
      text(output3$bifpoints[,1], output3$bifpoints[,9], output3$biftype, pos=2, offset=0.3)
    }
  }

  par(mfg=c(1,1), mar = c(0, 4, 0, 0.5))
  plot(1, 1, type="l", log="y", xaxt="n", xlim=c(0,0.0004), xlab="", xaxs="i",yaxs="i", yaxt="n", ylim=c(5.0E-6, 5.0E-4), ylab="")
  axis(1, at=(0:4)*0.0001, labels=c("", "", "", "", ""))
  axis(2, at=10^(-5:-4), labels=c(expression(10^{-5}), expression(10^{-4})), las=2)
  mtext("Predator density", 2, line=2.5, cex=1.3)

  if (stage > 2) {
    lines(output3$curvepoints[,1], output3$curvepoints[,3], type="l", col=rgb(0,.6,0), lwd=3)
    points(output3$bifpoints[,1], output3$bifpoints[,3], col="red", pch=8, lwd=2)
    text(output3$bifpoints[,1], output3$bifpoints[,3], output3$biftype, pos=2, offset=0.3)
  }

}

modelname="PNAS2002.h"

cat("\n\n\nStarting from the (known) trivial equilibrium with only resource\n\n")
cmd = paste0('output1 <- PSPMequi("', modelname, '", "EQ", c(1.0E-06, 1.0E-06), 0.5, c(1, 0, 4E-4), NULL, c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE)')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)
bifplot(1)

cat("\n\n\nStart from the detected branching point to compute the consumer-resource equilibrium\n\n")
cmd = paste0('output2 <- PSPMequi("', modelname, '", "EQ", output1$bifpoints[c(1, 2, 5)], 0.2, c(1, 0, 4E-4), NULL, c("envZE", "1", "envZE", "2"))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

cat('', '> output2$bifpoints', sep='\n'); print(output2$bifpoints); cat('', '> output2$biftypes', sep='\n'); print(output2$biftypes)
bifplot(2)

cat("\n\n\nStart from the detected branching point of the predator to compute the predator consumer-resource equilibrium\n\n")
cmd = paste0('output3 <- PSPMequi("', modelname, '", "EQ", output2$bifpoints[c(1, 2, 3, 7, 5)], -0.1, c(1, 0, 4E-4), NULL, NULL)')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

cat('', '> output3$bifpoints', sep='\n'); print(output3$bifpoints); cat('', '> output3$biftypes', sep='\n'); print(output3$biftypes)
bifplot(3)

cat("\n\n\nContinuation of the transcritical bifurcation curve of the structured population in 2 parameters (maximum resource density and background mortality)\n\n")
cmd = paste0('output4 <- PSPMequi("', modelname, '", "BP", c(output1$bifpoints[1:2], 0.01), 0.05, c(1, 0, 4E-4, 11, 0, 0.1), NULL, c("envZE", "1", "envZE", "2", "popBP", "0"))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

plot.new()
par(par.defaults)
par(cex = 1.0, cex.lab=1.3)
par(oma = c(0, 0, 0.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

plot(output4$curvepoints[,1], output4$curvepoints[,6], type="l", col=rgb(0,.6,0), lwd=3, xaxt="n", xlim=c(0,0.0004), xaxs="i", yaxs="i", yaxt="n", ylim=c(0.01, 0.09), xlab="Maximum resource density", ylab="")
axis(1, at=(0:4)*0.0001, labels=c("0.0", "0.0001", "0.0002", "0.0003", "0.0004"))
axis(2, at=c(1,3,5,7,9)*0.01, labels=TRUE, las=2)
mtext("Consumer mortality", 2, line=3.0, cex=1.3)

cat("\n\n\nContinuation of the transcritical bifurcation curve of the predator in 2 parameters (maximum resource density and background mortality)\n\n")
cmd = paste0('output5 <- PSPMequi("', modelname, '", "BPE", c(output2$bifpoints[c(1, 2, 5)], 0.01), -0.1, c(1, 0, 4E-4, 11, 0, 0.1), NULL, c("envZE", "2", "envBP", "1"))')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

lines(output5$curvepoints[,1], output5$curvepoints[,6], type="l", col=rgb(.6,0,0), lwd=3)

cat("\n\n\nContinuation of the saddle-node bifurcation in 2 parameters (maximum resource density and background mortality)\n\n")
cmd = paste0('output6 <- PSPMequi("', modelname, '", "LP", c(output3$bifpoints[1:5], 0.01), 0.05, c(1, 0, 4E-4, 11, 0, 0.1), NULL, NULL)')

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))

lines(output6$curvepoints[,1], output6$curvepoints[,6], type="l", col=rgb(0,0,0), lwd=3)

par(par.defaults)
PSPMclean('F')
setwd(oldwd)

