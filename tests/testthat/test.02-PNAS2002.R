modelname="PNAS2002.h"

# Start the clock!
ptm <- proc.time()

cmd = paste0('output1 <- PSPMequi("', modelname, '", "EQ", c(1.0E-06, 1.0E-06), 0.5, c(1, 0, 4E-4), NULL, c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE)')
eval(parse(text=cmd))

cat('', '> output1$bifpoints', sep='\n'); print(output1$bifpoints); cat('', '> output1$biftypes', sep='\n'); print(output1$biftypes)

cmd = paste0('output2 <- PSPMequi("', modelname, '", "EQ", output1$bifpoints[c(1, 2, 5)], 0.2, c(1, 0, 4E-4), NULL, c("envZE", "1", "envZE", "2"))')
eval(parse(text=cmd))

cat('', '> output2$bifpoints', sep='\n'); print(output2$bifpoints); cat('', '> output2$biftypes', sep='\n'); print(output2$biftypes)

cmd = paste0('output3 <- PSPMequi("', modelname, '", "EQ", output2$bifpoints[c(1, 2, 3, 7, 5)], -0.1, c(1, 0, 4E-4), NULL, NULL)')
eval(parse(text=cmd))

cat('', '> output3$bifpoints', sep='\n'); print(output3$bifpoints); cat('', '> output3$biftypes', sep='\n'); print(output3$biftypes)

cmd = paste0('output5 <- PSPMequi("', modelname, '", "BPE", c(output2$bifpoints[c(1, 2, 5)], 0.01), -0.1, c(1, 0, 4E-4, 11, 0, 0.1), NULL, c("envZE", "2", "envBP", "1"))')
eval(parse(text=cmd))

cmd = paste0('output6 <- PSPMequi("', modelname, '", "LP", c(output3$bifpoints[1:5], 0.01), 0.05, c(1, 0, 4E-4, 11, 0, 0.1), NULL, NULL)')
eval(parse(text=cmd))

# Do the numerical integration using the original EBT
initstate <- list(Environment = c(1.561276e-04, 1.270327e-04, 4.008016e-06), Pop00 = matrix(c(0.001, 0, 7.0, 1.0E-5, 300, 111), ncol = 3, byrow = TRUE))
PSPMecodyn(modelname, initstate, c(1, 1, 10, 500), options = c("report", "50"), force = TRUE)

# Do the numerical integration using the simplified EBT
CFLAGS="-DEBTMETHOD=0 "
initstate <- list(Environment = c(1.561276e-04, 1.270327e-04, 4.008016e-06), Pop00 = matrix(c(0.001, 0, 7.0, 1.0E-5, 300, 111), ncol = 3, byrow = TRUE))
PSPMecodyn(modelname, initstate, c(1, 1, 10, 500), options = c("report", "50"), force = TRUE)
if (exists("CFLAGS")) rm("CFLAGS")

# Stop the clock
tt <- proc.time() - ptm

print(tt)

PSPMclean("F")
