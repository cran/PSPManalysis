output <- PSPMdemo("Medfly", clean=TRUE, force=TRUE)

output <- PSPMdemo("Medfly", c(2, 11, 0.1, 11, 20))

# plot(output$curvepoints[,1], output$curvepoints[,2], type='l', lwd=2, xlab="Juvenile period", ylab="Population growth rate", font.lab=2)

PSPMclean("F")

