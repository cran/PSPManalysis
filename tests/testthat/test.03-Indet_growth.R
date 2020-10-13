modelname="Indet_growth.h"

cat('\n\n\nDetection of the ESS value of the bifurcation parameter of the structured population (ingestion exponent)\n\n');
cmd = paste0('output1 <- PSPMequi("', modelname, '", "EQ", c(1.0, 0.22, 0.0), -0.1, c(6, 0.5, 2.0), options=c("popEVO", "0", "report", "2"), clean = TRUE)')
eval(parse(text=cmd))

cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

cat('\n\n\nDetection of the ESS value of the bifurcation parameter (ingestion exponent) and the maintenance exponent of the structured population\n\n');
cmd = paste0('output1b <- PSPMequi("', modelname, '", "EQ", c(1.0, 0.22, 0.0), -0.1, c(6, 0.5, 2.0), options=c("popEVO", "0", "parEVO", "6", "popEVO", "0", "parEVO", "9", "report", "2"), clean = TRUE)')
eval(parse(text=cmd))

cat('', '> output1b$bifpoints', sep='\n'); print(output1b$bifpoints); cat('', '> output1b$biftypes', sep='\n'); print(output1b$biftypes)

cat("\n\n\nContinuation of the ESS value of the ingestion exponent of the structured population as a function of the first bifurcation parameter (maintenance exponent)\n\n")
cmd = c(paste0('output2a <- PSPMequi("', modelname, '", "ESS", c(1.0, output1$bifpoints[c(2, 3, 1)]), -0.1, c(9, 0.5, 2.0, 0, 6, 0.5, 2.0))'), paste0('output2b <- PSPMequi("', modelname, '", "ESS", c(1.0, output1$bifpoints[c(2, 3, 1)]), 0.1, c(9, 0.5, 2.0, 0, 6, 0.5, 2.0))')) 
eval(parse(text=cmd))

cat("\n\n\nConstruction of the PIP of the resident and mutant value of the first parameter of the structured population (ingestion exponent)\n\n")
cmd = c(paste0('output3 <- PSPMequi("', modelname, '", "PIP", c(output1$bifpoints[c(1, 2, 3, 1)]), 0.1, c(6, 0.5, 2.0, 6, 0.5, 2.0), options=c("popEVO", "0"))'), paste0('output4 <- PSPMequi("', modelname, '", "PIP", c(output1$bifpoints[c(1, 2, 3, 1)]), -0.1, c(6, 0.5, 2.0, 6, 0.5, 2.0), options=c("popEVO", "0"))'))
eval(parse(text=cmd))

cat("\n\n\nSimulating the simultaneous evolution of the ingestion and maintenance exponent over evolutionary time\n\n")
cmd = paste0('output4 <- PSPMevodyn("', modelname, '", c(0.22, 0.03554, 1.0, 1.0), c(0.05, 100), c(0, 6, 0.5, 1.5, 0, 9, 0.5, 1.5), options=c("report", "5"))')
eval(parse(text=cmd))


PSPMclean("F")

