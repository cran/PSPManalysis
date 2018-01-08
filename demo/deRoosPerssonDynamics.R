devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

# Do the numerical integration using the EBT: Trajectory 1
cat("\n\n\nSimulating the ecological dynamics of predators, consumers and resource from a computed equilibrium state\n\n")
str <- readline("Press any key to continue....\n")

output <- PSPMequi("PNAS2002", "EQ", c(3.0E-04, 1.561E-04, 1.270E-04, 4.008E-06, 2.761E-04), 0.1, c(1, 0, 1), options = c("single"), clean = TRUE, force = TRUE)
initstate <- csbread("PNAS2002-EQ-0000.csb", 1)
output1 <- PSPMecodyn("PNAS2002", initstate, c(1, 1, 10, 1000), options = c("report", "50"), clean = TRUE, force = TRUE)

plot.new()
par(par.defaults)
par(cex = 1.0, cex.lab=1.3)
par(oma = c(0, 0, 0.5, 2.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

plot(1, 1, type="l", xaxt="n", xlim=c(0,1000), xaxs="i", yaxs="i", log="y", yaxt="n", ylim=c(2.0E-6, 2.0E-3), xlab="Time", ylab="")
lines(output1$curvepoints[,1], output1$curvepoints[,7]+output1$curvepoints[,8], type="l", lwd=3, col=rgb(0,0,0.6))
lines(output1$curvepoints[,1], output1$curvepoints[,9], type="l", lwd=3, col=rgb(0.6,0,0))
axis(1, at=(0:5)*200, labels=TRUE)
axis(2, at=c(10,100,1000)*1.0E-6, labels=TRUE, las=2)
axis(2, at=c(1:9)*1.0E-6, labels=FALSE)
axis(2, at=c(1:9)*1.0E-5, labels=FALSE)
axis(2, at=c(1:9)*1.0E-4, labels=FALSE)
axis(2, at=c(1:9)*1.0E-3, labels=FALSE)
axis(4, at=c(10,100,1000)*1.0E-6, labels=TRUE, las=2)
axis(4, at=c(1:9)*1.0E-6, labels=FALSE)
axis(4, at=c(1:9)*1.0E-5, labels=FALSE)
axis(4, at=c(1:9)*1.0E-4, labels=FALSE)
axis(2, at=c(1:9)*1.0E-3, labels=FALSE)
mtext("Juvenile biomass", 2, line=3.0, cex=1.3, col=rgb(0,0,0.6))
mtext("Adult biomass", 4, line=3.0, cex=1.3, col=rgb(0.6,0,0))

# Do the numerical integration using the EBT: Trajectory 1
cat("\n\n\nSimulating the ecological dynamics of predators, consumers and resource from an initial state produced by PSPMind()\n\n")
str <- readline("Press any key to continue....\n")

initstate <- PSPMind("PNAS2002", c(1.561276e-04, 1.270327e-04, 4.008016e-06, 0.01), options = c("isort", "1"))
output2 <- PSPMecodyn("PNAS2002", initstate, c(1, 1, 10, 1000), options = c("report", "50"))

lines(output2$curvepoints[,1], output2$curvepoints[,7]+output2$curvepoints[,8], type="l", lty=2, lwd=3, col=rgb(0,0,0.6))
lines(output2$curvepoints[,1], output2$curvepoints[,9], type="l", lty=2, lwd=3, col=rgb(0.6,0,0))

# Do the numerical integration using the EBT: Trajectory 1
cat("\n\n\nSimulating the ecological dynamics of predators, consumers and resource from a manually constructed initial state\n\n")
str <- readline("Press any key to continue....\n")

initstate <- list(Environment = c(1.561276e-04, 1.270327e-04, 4.008016e-06), Pop00 = matrix(c(0.001, 0, 7.0, 1.0E-5, 300, 111), ncol = 3, byrow = TRUE))
output3 <- PSPMecodyn("PNAS2002", initstate, c(1, 1, 10, 1000), options = c("report", "50"))

lines(output3$curvepoints[,1], output3$curvepoints[,7]+output3$curvepoints[,8], type="l", lty=4, lwd=3, col=rgb(0,0,0.6))
lines(output3$curvepoints[,1], output3$curvepoints[,9], type="l", lty=4, lwd=3, col=rgb(0.6,0,0))

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
