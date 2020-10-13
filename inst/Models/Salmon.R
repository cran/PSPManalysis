PSPMdimensions <- c(PopulationNr = 1, IStateDimension = 2, 
                    LifeHistoryStages = 3, ImpactDimension = 5)

EnvironmentState <- c(X = "GENERALODE", P = "PERCAPITARATE")

DefaultParameters <- c(Rho = 0.01, Xmax = 0.5,
                       K = 1.0, Imax = 0.0025, Bmax = 0.002725,
                       L0 = 2.0, Ls = 20.0, Lm = 30.0,  Linf = 115.0,
                       Xi = 0.00051, Mu1 = 0.002, Mu2 = 0.006, 
                       D = 0.75, Phi = 0.001, Mup = 0.006)

StateAtBirth <- function(E, pars) {
  with(as.list(c(E, pars)),{
    c(Age = 0.0, Length = L0)
  })
}

LifeStageEndings <- function(lifestage, istate, birthstate, BirthStateNr, E, pars) {
  with(as.list(c(E, pars, istate)),{
    maturation = switch(lifestage, Length - Ls, Length - Lm, -1)
  })
}

LifeHistoryRates <- function(lifestage, istate, birthstate, BirthStateNr, E, pars) {
  with(as.list(c(E, pars, istate)),{
    list(
      development = c(1.0, 
                      switch(lifestage, Xi*(Linf*X/(K+X) - Length), 
                             Xi*(Linf - Length), Xi*(Linf - Length))),
      fecundity = switch(lifestage, 0, 0, Bmax*Length^2),
      mortality = switch(lifestage, Mu1, Mu2 + Phi*P*Length^(-D), 
                         Mu2 + Phi*P*Length^(-D)),
      impact = switch(lifestage, 
                      c(Imax*X/(K+X)*Length^2, 0, Length^3, 0, 0),
                      c(Imax*X/(K+X)*Length^2, Phi*Length^(3-D), 0, Length^3, 0),
                      c(Imax*X/(K+X)*Length^2, Phi*Length^(3-D), 0, 0, Length^3))
    )
  })
}

EnvEqui <- function(I, E, pars) {
  with(as.list(c(E, pars)),{
    c(Rho*(Xmax - X) - I[1], I[2] - Mup)
  })
}

