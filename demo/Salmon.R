devAskNewPage(ask = FALSE)
oldwd <- getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

modelfile <- "Salmon.h"

str <- readline(paste0("\nUse the C implementation of the Salmon model instead of its R implementation? [ Y / n ] :  "))
if ((str == "N") || (str == "n")) modelfile <- "Salmon.R"


###################################################################################
# Equilibrium predator, consumer and resource biomass densities as a function of  Xmax

# The function to plot the bifurcation plots in stages
Figure1a <- function(stage = 1) {
   devAskNewPage(ask = FALSE)
   # graphics.off()
   # plot.new()
   # Make the plot
   stableR   <- (EqR$curvepoints[,1] <= EqR$bifpoints[1,1])
   if (stage > 1) stableCR  <- (EqCR$curvepoints[,1] <= EqCR$bifpoints[1,1])
   if (stage > 2) stablePCR <- (EqPCR$curvepoints[,2] >= EqPCR$bifpoints[1,2])
   
   layout(matrix((1:3), nrow = 3, ncol = 1), heights = c(0.85, 0.85, 1))
   par(tcl = 0.5)
   par(mar = c(0, 10, 2, 10))
   cexbase <- dev.size()[2]/16.6
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n", 
        xlim=c(0.0, 10), ylim = c(0, 200), xlab="", ylab="")
   
   lines(c(0,EqR$bifpoints[,1]), c(0, 0), lwd = 2)
   points(EqR$bifpoints[,1], 0, col="red", pch=8, lwd=2, cex = 2)
   text(EqR$bifpoints[,1], 0, EqR$biftype, pos=3, offset=0.7, cex = 1.5)
   if (stage > 1) {
      lines(c(EqR$bifpoints[,1], EqCR$bifpoints[,1]), c(0, 0), lwd = 4)
      points(EqCR$bifpoints[,1], EqCR$bifpoints[,3], col="red", pch=8, lwd=2, cex = 2)
      text(EqCR$bifpoints[,1], EqCR$bifpoints[,3], EqCR$biftype, pos=4, offset=0.7, cex = 1.5)
   }
   if (stage > 2) {
      lines(EqPCR$curvepoints[stablePCR,1], EqPCR$curvepoints[stablePCR,3], lwd = 6)
      lines(EqPCR$curvepoints[!stablePCR,1], EqPCR$curvepoints[!stablePCR,3], lwd = 6, lty = "dashed")
      points(EqPCR$bifpoints[,1], EqPCR$bifpoints[,3], col="red", pch=8, lwd=2, cex = 2)
      text(EqPCR$bifpoints[,1], EqPCR$bifpoints[,3], EqPCR$biftype, pos=2, offset=0.7, cex = 1.5)
   }
   axis(1, at = (0:6)*2, labels=F)
   axis(2, at = (0:6)*40, labels=T, las=2, cex.axis = 1.6)
   mtext(expression("Predators" ~ ("#" %.% m^{-3})), 2, line = 6, cex = 1.8*cexbase)
   
   par(mar = c(0, 10, 0, 10))
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n", xlim=c(0.0, 10), ylim = c(0, 1300), xlab="", ylab="")
   lines(c(0, EqR$bifpoints[,1]), c(0, 0), lwd = 2, col = rgb(0, 0, 0.6))
   lines(c(0, EqR$bifpoints[,1]), c(0, 0), lwd = 2, col = rgb(0.6, 0, 0))
   points(EqR$bifpoints[,1], 0, col="red", pch=8, lwd=2, cex = 2)
   text(EqR$bifpoints[,1], 0, EqR$biftype, pos=3, offset=0.7, cex = 1.5)
   
   lraxisratio <- 10.0
   if (stage > 1) {
      lines(EqCR$curvepoints[stableCR,1], EqCR$curvepoints[stableCR,7], lwd = 4, col = rgb(0, 0, 0.6))
      points(EqCR$bifpoints[,1], EqCR$bifpoints[,7], col="red", pch=8, lwd=2, cex = 2)
      text(EqCR$bifpoints[,1], EqCR$bifpoints[,7], EqCR$biftype, pos=4, offset=0.7, cex = 1.5)

      lines(EqCR$curvepoints[stableCR,1], lraxisratio*(EqCR$curvepoints[stableCR,8]+EqCR$curvepoints[stableCR,9]), lwd = 4, col = rgb(0.6, 0, 0))
      points(EqCR$bifpoints[,1], lraxisratio*(EqCR$bifpoints[,8]+EqCR$bifpoints[,9]), col="red", pch=8, lwd=2, cex = 2)
      text(EqCR$bifpoints[,1], lraxisratio*(EqCR$bifpoints[,8]+EqCR$bifpoints[,9]), EqCR$biftype, pos=4, offset=0.7, cex = 1.5)
      
      # indx <- which.min(abs(EqCR$curvepoints[,1] - 4.0))
      # points(4.0, EqCR$curvepoints[indx,7], col=gray(0.6), pch=16, cex = 3)
      # points(4.0, lraxisratio*(EqCR$curvepoints[indx,8]+EqCR$curvepoints[indx,9]), col=gray(0.6), pch=16, cex = 3)
   }
   if (stage > 2) {
      lines(EqPCR$curvepoints[stablePCR,1], EqPCR$curvepoints[stablePCR,7], lwd = 6, col = rgb(0, 0, 0.6))
      lines(EqPCR$curvepoints[!stablePCR,1], EqPCR$curvepoints[!stablePCR,7], lwd = 6, col = rgb(0, 0, 0.6), lty = "dashed")
      points(EqPCR$bifpoints[,1], EqPCR$bifpoints[,7], col="red", pch=8, lwd=2, cex = 2)
      text(EqPCR$bifpoints[,1], EqPCR$bifpoints[,7], EqPCR$biftype, pos=1, offset=1.0, cex = 1.5)

      lines(EqPCR$curvepoints[stablePCR,1], lraxisratio*(EqPCR$curvepoints[stablePCR,8]+EqPCR$curvepoints[stablePCR,9]), lwd = 6, col = rgb(0.6, 0, 0))
      lines(EqPCR$curvepoints[!stablePCR,1], lraxisratio*(EqPCR$curvepoints[!stablePCR,8]+EqPCR$curvepoints[!stablePCR,9]), lwd = 6, col = rgb(0.6, 0, 0), lty = "dashed")
      points(EqPCR$bifpoints[,1], lraxisratio*(EqPCR$bifpoints[,8]+EqPCR$bifpoints[,9]), col="red", pch=8, lwd=2, cex = 2)
      text(EqPCR$bifpoints[,1], lraxisratio*(EqPCR$bifpoints[,8]+EqPCR$bifpoints[,9]), EqPCR$biftype, pos=3, offset=1.0, cex = 1.5)
      
      # y <- EqPCR$curvepoints[stablePCR, c(1, 7, 8, 9)]
      # indx <- which.min(abs(y[,1] - 4.0))
      # points(4.0, y[indx,2], col=gray(0.0), pch=16, cex = 3)
      # points(4.0, lraxisratio*(y[indx,3]+y[indx,4]), col=gray(0.0), pch=16, cex = 3)
   }

   axis(1, at = (0:6)*2, labels=F)
   axis(2, at = (0:6)*200, labels=T, las=2, cex.axis = 1.6)
   rticks <- (0:6)*20
   axis(4, at = rticks*lraxisratio, labels=rticks, las=2, cex.axis = 1.6)
   mtext(expression("Small consumers" ~ (g %.% m^{-3})), 2, line = 6, cex = 1.8*cexbase)
   mtext(expression("Large consumers" ~ (g %.% m^{-3})), 4, line = 6, cex = 1.8*cexbase)
   
   par(mar = c(7, 10, 0, 10))
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n", 
        xlim=c(0.0, 10), ylim = c(0, 9), xlab="", ylab="")
   lines(c(0, EqR$curvepoints[stableR,1]), c(0,EqR$curvepoints[stableR,2]), lwd = 2, col = rgb(0, 0.6, 0))
   points(EqR$bifpoints[,1], EqR$bifpoints[,2], col="red", pch=8, lwd=2, cex = 2)
   text(EqR$bifpoints[,1], EqR$bifpoints[,2], EqR$biftype, pos=3, offset=0.7, cex = 1.5)
   
   if (stage > 1) {
      lines(EqCR$curvepoints[stableCR,1], EqCR$curvepoints[stableCR,2], lwd = 4, col = rgb(0, 0.6, 0))
      points(EqCR$bifpoints[,1], EqCR$bifpoints[,2], col="red", pch=8, lwd=2, cex = 2)
      text(EqCR$bifpoints[,1], EqCR$bifpoints[,2], EqCR$biftype, pos=4, offset=0.7, cex = 1.5)
   }
   if (stage > 2) {
      lines(EqPCR$curvepoints[stablePCR,1], EqPCR$curvepoints[stablePCR,2], lwd = 6, col = rgb(0, 0.6, 0))
      lines(EqPCR$curvepoints[!stablePCR,1], EqPCR$curvepoints[!stablePCR,2], lwd = 6, col = rgb(0, 0.6, 0), lty = "dashed")
      points(EqPCR$bifpoints[,1], EqPCR$bifpoints[,2], col="red", pch=8, lwd=2, cex = 2)
      text(EqPCR$bifpoints[,1], EqPCR$bifpoints[,2], EqPCR$biftype, pos=3, offset=1.0, cex = 1.5)
   }
   axis(1, at = (0:6)*2, labels=T, cex.axis = 1.6)
   axis(2, at = (0:6)*2, labels=T, las=2, cex.axis = 1.6)
   mtext(expression("Maximum resource density" ~ (g %.% m^{-3})), 1, line = 5, cex = 1.8*cexbase)
   mtext(expression("Resource" ~ (g %.% m^{-3})), 2, line = 6, cex = 1.8*cexbase)
}

cat("\n#############################################################################################\n\n")
cat("Figure 1a: Complete bifurcation graph with equilibrium predator, consumer and resource\n")
cat("           biomass densities as a function of the maximum resource density\n")

########################### Only resources
cat("\n")
cat("Computation step 1: Starting from the (known) trivial equilibrium with only resource,\n")
cat("                    compute a curve of resource-only steady states\n\n")

cmd <- paste0('EqR <- PSPMequi(modelname = "', modelfile, '", biftype = "EQ", startpoint = c(0.1, 0.1), 
                stepsize =  0.5, parbnds = c(1, 0.01, 10), options = c("popZE", "0", "envZE", "1"), clean = TRUE)')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure1a(1)

########################### Consumers and resources
cat("#############################################################################################\n\n")
cat("Computation step 2: Starting from the detected bifurcation (branching) point,\n")
cat("                    compute a curve of consumer-resource steady states\n\n")

cmd <- paste0('EqCR <- PSPMequi(modelname = "', modelfile, '", biftype = "EQ", startpoint = c(EqR$bifpoints[1,1:2], 0), 
                 stepsize =  0.5, parbnds = c(1, 0.01, 10), options = c("envZE", "1"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure1a(2)

########################### Predators, consumers and resources
cat("#############################################################################################\n\n")
cat("Computation step 3: Starting from the detected bifurcation (branching) point,\n")
cat("                    compute a curve of predator-consumer-resource steady states\n\n")

cmd <- paste0('EqPCR <- PSPMequi(modelname = "', modelfile, '", biftype = "EQ", startpoint = EqCR$bifpoints[1,1:4], 
                  stepsize = -0.5, parbnds = c(1, 0.01, 10))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure1a(3)

###################################################################################
# Population state plots at Xmax = 4.0 in predator-consumer-resource and consumer-resource equilibrium

Figure1b <- function(stage = 1) {
   par(par.defaults)

   layout(matrix((1:2), nrow = 2, ncol = 1))
   par(tcl = 0.5)
   par(mar = c(6, 10, 3, 10))
   cexbase <- dev.size()[2]/16.6
   
   binwidth <- 5.0
   xmax <- 65
   ymin <- 1.0E-6
   ymax <- 1.0
   
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n", xlab="", ylab="", log = "y", xlim = c(0,xmax), ylim = c(ymin,ymax))
   axis(1, at = (0:6)*10, labels=T, cex.axis = 1.6*cexbase)
   axis(2, at = 10^seq(-6, 0, 2), labels=c(expression(10^{-6}), expression(10^{-4}), expression(10^{-2}), 1.0), las=2, cex.axis = 1.6*cexbase)
   mtext("Body length (cm)", 1, line = 3, cex = 1.8*cexbase)
   mtext(expression("Consumer density" ~ ("#" %.% m^{-3})), 2, line = 5, cex = 1.8*cexbase)
   y <- psCR$Pop00[,c(3,1)]
   binindex <- ceiling(y[,1] / binwidth)
   binyvals <- unlist(lapply((1:(xmax/binwidth)), function(indx) {sum(y[binindex == indx,2])} ))
   x <- ((1:length(binyvals)) - 0.25 ) * binwidth
   rect(xleft=x-binwidth/4.0, ybottom=ymin, xright=x+binwidth/4.0, ytop=pmax(binyvals, ymin), col=gray(0.6))
   
   if (stage > 1) {
      y <- psPCR[[2]]$Pop00[,c(3,1)]
      binindex <- ceiling(y[,1] / binwidth)
      binyvals <- unlist(lapply((1:(xmax/binwidth)), function(indx) {sum(y[binindex == indx,2])} ))
      x <- ((1:length(binyvals)) - 0.75 ) * binwidth
      rect(xleft=x-binwidth/4.0, ybottom=ymin, xright=x+binwidth/4.0, ytop=pmax(binyvals, ymin), col=gray(0.0))
   }

   par(mar = c(8, 10, 1, 10))
   plot(psCR$Pop00[,2], psCR$Pop00[,3], type="l", xaxt="n", yaxt="n", xlab="", ylab="", xlim = c(0,4000), ylim = c(0,100), lwd = 3, col = gray(0.6))
   axis(1, at = (0:4)*1000, labels=T, cex.axis = 1.6*cexbase)
   axis(2, at = (0:5)*20, labels=T, cex.axis = 1.6*cexbase, las = 2)
   mtext("Age (days)", 1, line = 3, cex = 1.8*cexbase)
   mtext("Body length (cm)", 2, line = 5, cex = 1.8*cexbase)
   if (stage > 1) {
      lines(psPCR[[2]]$Pop00[,2], psPCR[[2]]$Pop00[,3], lwd = 3, col = gray(0.0))
   }
}
   
cat("\n#############################################################################################\n\n")
cat("Figure 1b: Plot the length distribution and the length-age relation of the consumer population in the stable consumer-resource\n")
cat("           and the stable predator-consumer-resource equilibrium at maximum resource density equal to 4.0\n")

cat("\n")
cat("Step 1: Extract the consumer population state in the stable consumer-resource equilibrium\n")
cat("        at maximum resource density equal to 4.0 and plot the population length distribution\n")
cat("        and the length-age relationship\n\n")

cmd <- paste0('psCR <- csbread("Salmon-EQ-0001.csb", "State-4.0")')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure1b(1)

cat("\n")
cat("Step 2: Extract the consumer population state in the stable predator-consumer-resource equilibrium\n")
cat("        at maximum resource density equal to 4.0 and plot the population length distribution\n")
cat("        and the length-age relationship\n\n")

cmd <- paste0('psPCR <- csbread("Salmon-EQ-0002.csb", "State-4.0")')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure1b(2)

###################################################################################
# Two parameter plot of LP and BPE as a function of Xmax and Mup

# The function to plot the bifurcation plots in stages
Figure2 <- function(stage = 1) {
   devAskNewPage(ask = FALSE)
   # graphics.off()
   # plot.new()
   # Make the plot
   par(par.defaults)
   par(mar = c(6, 10, 1, 10))
   plot(NULL, NULL, type = "n", xaxt="n", yaxt="n", xlim=c(0.0, 8), ylim = c(0, 0.012), xlab="", ylab="")
   lines(BPEm$curvepoints[,1], BPEm$curvepoints[,5], lwd = 3, col = rgb(0, 0, 0.6))
   if (stage > 1) lines(BPEp$curvepoints[,1], BPEp$curvepoints[,5], lwd = 3, col = rgb(0, 0, 0.6))
   if (stage > 2) lines(LPm$curvepoints[,1], LPm$curvepoints[,5], lwd = 3, col = rgb(0.6, 0, 0))
   if (stage > 3) lines(LPp$curvepoints[,1], LPp$curvepoints[,5], lwd = 3, col = rgb(0.6, 0, 0))
   axis(1, at = (0:6)*2, labels=T, cex.axis = 1.6)
   axis(2, at = (0:6)*0.002, labels=T, las=2, cex.axis = 1.6)
   mtext(expression("Maximum resource density" ~ (g %.% m^{-3})), 1, line = 4, cex = 1.8)
   mtext(expression("Predator mortality" ~ (day^{-1})), 2, line = 6, cex = 1.8)
   legend("topleft", c("BPE #1", "LP"),  col = c(rgb(0, 0, 0.6), rgb(0.6, 0, 0)), lwd = 3, cex = 1.5)
}

cat("\n#############################################################################################\n\n")
cat("Figure 2: Location of the predator invasion boundary (BPE #1 in Figure 1) dependent on\n")
cat("          the maximum resource density and predator mortality\n")

cat("\n")
cat("Computation step 1: Starting from the bifurcation point BPE #1 detected in Figure 1,\n")
cat("                    compute the predator invasion boundary for decreasing maximum resource densities\n\n")

cmd <- paste0('BPEm <- PSPMequi(modelname = "', modelfile, '", biftype = "BPE", startpoint = c(EqCR$bifpoints[1,c(1:2,4)], 0.006), 
                 stepsize = -0.5, parbnds = c(1, 0.0, 10, 14, 0, 0.05), options = c("envBP", "1"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure2(1)

cat("#############################################################################################\n\n")
cat("Computation step 2: Starting from the bifurcation point BPE #1 detected in Figure 1,\n")
cat("                    compute the predator invasion boundary for increasing maximum resource densities\n\n")

cmd <- paste0('BPEp <- PSPMequi(modelname = "', modelfile, '", biftype = "BPE", startpoint = c(EqCR$bifpoints[1,c(1:2,4)], 0.006), 
                 stepsize = 0.5, parbnds = c(1, 0.0, 10, 14, 0, 0.05), options = c("envBP", "1"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure2(2)

cat("#############################################################################################\n\n")
cat("Computation step 3: Starting from the bifurcation point LP detected in Figure 1,\n")
cat("                    compute the predator persistence boundary for decreasing maximum resource densities\n\n")

cmd <- paste0('LPm <- PSPMequi(modelname = "', modelfile, '", biftype = "LP", startpoint = c(EqPCR$bifpoints[1,1:4], 0.006), 
                stepsize = -0.5, parbnds = c(1, 0.0, 10, 14, 0, 0.05))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure2(3)

cat("#############################################################################################\n\n")
cat("Computation step 4: Starting from the bifurcation point LP detected in Figure 1,\n")
cat("                    compute the predator persistence boundary for increasing maximum resource densities\n\n")

cmd <- paste0('LPp <- PSPMequi(modelname = "', modelfile, '", biftype = "LP", startpoint = c(EqPCR$bifpoints[1,1:4], 0.006), 
                stepsize = 0.5, parbnds = c(1, 0.0, 10, 14, 0, 0.05))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure2(4)

###################################################################################
# Equilibrium predator, consumer and resource biomass densities as a function of Ls

# The function to plot the bifurcation plots in stages
Figure3 <- function(stage = 1) {
   devAskNewPage(ask = FALSE)
   # graphics.off()
   # plot.new()
   # Make the plot
   # Make the plot
   stableCR <-  (EvoCR$curvepoints[,1]  >= EvoCR$bifpoints[1,1])
   if (stage > 1) stablePCR <- (EvoPCR$curvepoints[,2] >= EvoPCR$bifpoints[1,2])
   
   layout(matrix((1:3), nrow = 3, ncol = 1), heights = c(0.85, 0.85, 1))
   par(tcl = 0.5)
   par(mar = c(0, 10, 2, 10))
   cexbase <- dev.size()[2]/16.6
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n", xlim=c(5.0, 9), ylim = c(0, 45), xlab="", ylab="")

   lines(c(0,EvoCR$bifpoints[,1]), c(0, 0), lwd = 3, lty = "dashed")
   lines(c(EvoCR$bifpoints[,1], 25), c(0, 0), lwd = 3)
   points(EvoCR$bifpoints[,1], EvoCR$bifpoints[,3], col="red", pch=8, lwd=2, cex = 2)
   text(EvoCR$bifpoints[,1], EvoCR$bifpoints[,3], EvoCR$biftype, pos=3, offset=1.5, cex = 1.5)
   
   if (stage > 1) {
      lines(EvoPCR$curvepoints[stablePCR,1], EvoPCR$curvepoints[stablePCR,3], lwd = 3)
      lines(EvoPCR$curvepoints[!stablePCR,1], EvoPCR$curvepoints[!stablePCR,3], lwd = 3, lty = "dashed")
      
      points(EvoPCR$bifpoints[1,1], EvoPCR$bifpoints[1,3], col="red", pch=8, lwd=2, cex = 2)
      text(EvoPCR$bifpoints[1,1], EvoPCR$bifpoints[1,3], EvoPCR$biftype[1], pos=2, offset=0.9, cex = 1.5)
      points(EvoPCR$bifpoints[2,1], EvoPCR$bifpoints[2,3], col="red", pch=8, lwd=2, cex = 2)
      text(EvoPCR$bifpoints[2,1], EvoPCR$bifpoints[2,3], EvoPCR$biftype[2], pos=3, offset=0.9, cex = 1.5)
   }
   axis(1, at = (0:10)*1, labels=F)
   axis(2, at = (0:6)*10, labels=T, las=2, cex.axis = 1.6)
   mtext(expression("Predators" ~ ("#" %.% m^{-3})), 2, line = 6, cex = 1.8*cexbase)
   
   par(mar = c(0, 10, 0, 10))
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n", xlim=c(5.0, 9), ylim = c(0, 75), xlab="", ylab="")
   lines(EvoCR$curvepoints[stableCR,1], EvoCR$curvepoints[stableCR,7], lwd = 3, col = rgb(0, 0, 0.6))
   points(EvoCR$bifpoints[,1], EvoCR$bifpoints[,7], col="red", pch=8, lwd=2, cex = 2)
   text(EvoCR$bifpoints[,1], EvoCR$bifpoints[,7], EvoCR$biftype, pos=2, offset=1.0, cex = 1.5)
   
   lraxisratio <- 0.7
   if (stage > 1) {
      lines(EvoPCR$curvepoints[stablePCR,1], EvoPCR$curvepoints[stablePCR,7], lwd = 3, col = rgb(0, 0, 0.6))
      lines(EvoPCR$curvepoints[!stablePCR,1], EvoPCR$curvepoints[!stablePCR,7], lwd = 3, col = rgb(0, 0, 0.6), lty = "dashed")
      points(EvoPCR$bifpoints[1,1], EvoPCR$bifpoints[1,7], col="red", pch=8, lwd=2, cex = 2)
      text(EvoPCR$bifpoints[1,1], EvoPCR$bifpoints[1,7], EvoPCR$biftype[1], pos=4, offset=1.3, cex = 1.5)
      points(EvoPCR$bifpoints[2,1], EvoPCR$bifpoints[2,7], col="red", pch=8, lwd=2, cex = 2)
      text(EvoPCR$bifpoints[2,1], EvoPCR$bifpoints[2,7], EvoPCR$biftype[2], pos=1, offset=1.3, cex = 1.5)
   }
   
   lines(EvoCR$curvepoints[stableCR,1], lraxisratio*(EvoCR$curvepoints[stableCR,8]+EvoCR$curvepoints[stableCR,9]), lwd = 3, col = rgb(0.6, 0, 0))
   points(EvoCR$bifpoints[,1], lraxisratio*(EvoCR$bifpoints[,8]+EvoCR$bifpoints[,9]), col="red", pch=8, lwd=2, cex = 2)
   text(EvoCR$bifpoints[,1], lraxisratio*(EvoCR$bifpoints[,8]+EvoCR$bifpoints[,9]), EvoCR$biftype, pos=2, offset=1.0, cex = 1.5)

   if (stage > 1) {
      lines(EvoPCR$curvepoints[stablePCR,1], lraxisratio*(EvoPCR$curvepoints[stablePCR,8]+EvoPCR$curvepoints[stablePCR,9]), lwd = 3, col = rgb(0.6, 0, 0))
      lines(EvoPCR$curvepoints[!stablePCR,1], lraxisratio*(EvoPCR$curvepoints[!stablePCR,8]+EvoPCR$curvepoints[!stablePCR,9]), lwd = 3, col = rgb(0.6, 0, 0), lty = "dashed")
      points(EvoPCR$bifpoints[,1], lraxisratio*(EvoPCR$bifpoints[,8]+EvoPCR$bifpoints[,9]), col="red", pch=8, lwd=2, cex = 2)
      text(EvoPCR$bifpoints[1,1], lraxisratio*(EvoPCR$bifpoints[1,8]+EvoPCR$bifpoints[1,9]), EvoPCR$biftype[1], pos=1, offset=1.0, cex = 1.5)
      text(EvoPCR$bifpoints[2,1], lraxisratio*(EvoPCR$bifpoints[2,8]+EvoPCR$bifpoints[2,9]), EvoPCR$biftype[2], pos=3, offset=1.0, cex = 1.5)
   }

   axis(1, at = (0:10)*1, labels=F)
   axis(2, at = (0:6)*20, labels=T, las=2, cex.axis = 1.6)
   rticks <- (0:6)*20
   axis(4, at = rticks*lraxisratio, labels=rticks, las=2, cex.axis = 1.6)
   mtext(expression("Small consumers" ~ (g %.% m^{-3})), 2, line = 6, cex = 1.8*cexbase)
   mtext(expression("Large consumers" ~ (g %.% m^{-3})), 4, line = 6, cex = 1.8*cexbase)

   par(mar = c(10, 10, 0, 10))
   plot(NULL, NULL, type="n", xaxt="n", yaxt="n",  xlim=c(5.0, 9), ylim = c(-1.2, 0.3), xlab="", ylab="")
   lines(EvoCR$curvepoints[stableCR,1], EvoCR$curvepoints[stableCR,12], lwd = 3, col = rgb(0, 0.6, 0))
   lines(par("usr")[1:2], c(0,0), lwd = 1, lty = "dashed")
   points(EvoCR$bifpoints[,1], EvoCR$bifpoints[,12], col="red", pch=8, lwd=2, cex = 2)
   text(EvoCR$bifpoints[,1], EvoCR$bifpoints[,12], EvoCR$biftype, pos=2, offset=1.0, cex = 1.5)

   if (stage > 1) {
      lines(EvoPCR$curvepoints[stablePCR,1], EvoPCR$curvepoints[stablePCR,12], lwd = 3, col = rgb(0, 0.6, 0))
      lines(EvoPCR$curvepoints[!stablePCR,1], EvoPCR$curvepoints[!stablePCR,12], lwd = 3, col = rgb(0, 0.6, 0), lty = "dashed")
      points(EvoPCR$bifpoints[,1], EvoPCR$bifpoints[,12], col="red", pch=8, lwd=2, cex = 2)
      text(EvoPCR$bifpoints[1,1], EvoPCR$bifpoints[1,12], EvoPCR$biftype[1], pos=2, offset=1.0, cex = 1.5)
      text(EvoPCR$bifpoints[2,1], EvoPCR$bifpoints[2,12], EvoPCR$biftype[2], pos=1, offset=1.0, cex = 1.5)
   }
   axis(1, at = (0:10)*1, labels=T, cex.axis = 1.6)
   axis(2, at = -1.2 + (0:3)*0.4, labels=c("-1.2", "-0.8", "-0.4", "0"), las=2, cex.axis = 1.6)
   mtext("Body size at habitat shift (cm)", 1, line = 5, cex = 1.8*cexbase)
   mtext("Selection gradient", 2, line = 7, cex = 1.8*cexbase)
   mtext(expression("(d" * italic(R)[0] * "/d" * italic(l)[s] * ")" ~ "(-)"), 2, line = 4, cex = 1.7*cexbase)
}

cat("\n#############################################################################################\n\n")
cat("Figure 3: Complete bifurcation graph with equilibrium predator, consumer and resource\n")
cat("          biomass densities as a function of body size at habitat shift\n")

########################### Consumers and resources
cat("\n")
cat("Computation step 1: Compute the curve of consumer-resource steady states\n")
cat("                    for decreasing values of body size at habitat shift\n\n")

cmd <- paste0('EvoCR <- PSPMequi(modelname = "', modelfile, '", biftype = "EQ", startpoint = c(25, 3.01798283E-01, 1.05872962E-04),
                  stepsize = -0.1, parbnds = c(6, 5.0, 25), options = c("envZE", "1", "popEVO", "0"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure3(1)

########################### Predators, consumers and resources
cat("#############################################################################################\n\n")
cat("Computation step 2: Starting from the detected bifurcation (branching) point,\n")
cat("                    compute a curve of predator-consumer-resource steady states\n\n")

cmd <- paste0('EvoPCR <- PSPMequi(modelname = "', modelfile, '", biftype = "EQ", startpoint = EvoCR$bifpoints[1,1:4], 
                   stepsize = 0.2, parbnds = c(6, 5.0, 25), options = c("popEVO", "0"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure3(2)

###################################################################################
# PIP as a function of resident and mutant Ls value

# The function to plot the bifurcation plots in stages
Figure4a <- function(stage = 1) {
   devAskNewPage(ask = FALSE)
   # graphics.off()
   # plot.new()
   # Make the plot
   par(par.defaults)
   par(mar = c(6, 10, 1, 10))
   plot(NULL, NULL, type = "n", xaxt="n", yaxt="n",  xlim=c(4, 9), ylim = c(4, 9), xlab="", ylab="")
   lines(PIPp$curvepoints[,1], PIPp$curvepoints[,5], lwd = 3, col = rgb(0, 0, 0.6))
   if (stage > 1) lines(PIPm$curvepoints[,1], PIPm$curvepoints[,5], lwd = 3, col = rgb(0, 0, 0.6))
   lines(par("usr")[1:2], par("usr")[3:4], lwd = 1, col = "black", lty = "solid")
   axis(1, at = (0:10)*1, labels=T, cex.axis = 1.6)
   axis(2, at = (0:10)*1, labels=T, las=2, cex.axis = 1.6)
   mtext("Resident body size at habitat shift (cm)", 1, line = 4, cex = 1.8)
   mtext("Mutant body size at habitat shift value (cm)", 2, line = 5, cex = 1.8)
}

cat("\n#############################################################################################\n\n")
cat("Figure 4a: Pairwise invasibility plot as function the resident and mutant value of\n")
cat("           the body size at habitat shift\n")

cat("\n")
cat("Computation step 1: Starting from the evolutionary steady state (ESS) detected in Figure 3, compute the isocline\n")
cat("                    with zero mutant fitness for increasing resident body sizes at habitat shift\n\n")

cmd <- paste0('PIPp <- PSPMequi(modelname = "', modelfile, '", biftype = "PIP", startpoint = EvoPCR$bifpoints[2,c(1:4,1)], 
                 stepsize = 0.1, parbnds = c(6, 2, 15, 6, 2, 15), options = c("popEVO", "0"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure4a(1)

cat("#############################################################################################\n\n")
cat("Computation step 2: Starting from the evolutionary steady state (ESS) detected in Figure 3, compute the isocline\n")
cat("                    with zero mutant fitness for decreasing resident body sizes at habitat shift\n\n")

cmd <- paste0('PIPm <- PSPMequi(modelname = "', modelfile, '", biftype = "PIP", startpoint = EvoPCR$bifpoints[2,c(1:4,1)], 
                 stepsize = -0.1, parbnds = c(6, 2, 15, 6, 2, 15), options = c("popEVO", "0"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

Figure4a(2)

###################################################################################
# Evolutionary dynamics of Ls

cat("\n#############################################################################################\n\n")
cat("Figure 4b: Dynamics of the value of the length at habitat shift over evolutionary time\n")
cat("           as predicted by the canonical equation of Adaptive Dynamics\n\n")

cmd <- paste0('TSevo <- PSPMevodyn(modelname = "', modelfile, '", startpoint = c(0.2566152, 35.55643, 0.0006141598, 5.019757), 
                    curvepars = c(10,100000), evopars = c(0, 6, 5, 10), options = c("report", "100"))')
cat(cmd)
str <- readline("\n\nPress any key to continue....\n")
eval(parse(text=cmd))

# Plot only the dynamics of Ls over evolutionary time
par(par.defaults)
par(mar = c(6, 10, 1, 10))
plot(NULL, NULL, type="n", xaxt="n", yaxt="n", 
     xlim=c(0, 8.0E+4), ylim = c(4.8, 6.3), xlab="", ylab="")
lines(TSevo$curvepoints[,1], TSevo$curvepoints[,5], lwd = 3, col = rgb(0, 0.6, 0))
lines(par("usr")[1:2], EvoPCR$bifpoints[2,c(1,1)], lwd = 1, lty = "dashed")
axis(1, at = (0:4)*2.0E+4, 
     labels=c("0", "20000", "40000", "60000", "80000"), cex.axis = 1.6)
axis(2, at = 5 + (0:10)*0.4, labels=T, las=2, cex.axis = 1.6)
mtext("Evolutionary time (scaled units)", 1, line = 4, cex = 1.8)
mtext("Body size at habitat shift (cm)", 2, line = 5, cex = 1.8)

###################################################################################
# Clean up
PSPMclean("F")

setwd(oldwd)
par(par.defaults)
rm(list = c("modelfile", "oldwd", "str", "cmd", "BPEm", "BPEp", "EqR", "EqCR", "EqPCR", "EvoCR", "EvoPCR", "LPm", "LPp", "PIPm", "PIPp", "TSevo", "par.defaults",  "psCR", "psPCR", "Figure1a", "Figure1b", "Figure2", "Figure3", "Figure4a"))

