devAskNewPage(ask = FALSE)
oldwd = getwd()
setwd(system.file("Models", package="PSPManalysis"))
if (!exists("par.defaults")) par.defaults <- par(no.readonly = T)

cat('\n\n\nContinuation of the equilibrium in the Kooijman DEB model as a function of a toxic maintenance effect (Martins et al. 2013, Am. Nat.)\n\n');
cmd = 'output1 <- PSPMequi("MartinsDEB", "EQ", c(1.31634938, 10.0, 0.390274776, -6.89E-07), -0.5, c(18, 0.0, 1.5), clean=TRUE, force=TRUE)'

str = readline(paste0("\n> ", cmd, "\n\nPress any key to continue....\n"))
eval(parse(text=cmd))
cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)
  
plot(output1$curvepoints[,1], output1$curvepoints[,10], type='l', col=rgb(0,0,0.6), lwd=3, xlim=c(0.0,1.4), ylim=c(0.0,10000), xlab="Toxic effect size", ylab="Juvenile & adult biomass", font.lab=2)
lines(output1$curvepoints[,1], output1$curvepoints[,11], type='l', col=rgb(0.6,0,0), lwd=3)

par(par.defaults)
PSPMclean('F')
setwd(oldwd)
