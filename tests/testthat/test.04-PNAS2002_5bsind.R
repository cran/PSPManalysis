pars = c(0.1, 8.48186E-05, 7.0, 27.0, 110.0, 300.0, 9.0E-6, 1.0E-4, 1.5E-5, 0.006, 0.003, 0.01, 5000.0, 0.1, 0.5, 0.01);

output <- PSPMind("PNAS2002_5bs", c(1.30341E-05, 3.84655E-05, 4.00802E-06), pars, options = c("isort", "1"), clean=TRUE, force=TRUE)

print(output)

PSPMclean("F")

