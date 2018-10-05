rm(list=ls())
library(BIOMASS)
library(profvis)
library(microbenchmark)

data("KarnatakaForest")
data(NouraguesHD)

myrtruncnorm <- function(n,lower = -1, upper = 1,mean=0,sd=1) {
  qnorm(runif(n,pnorm(lower,mean=mean,sd=sd),pnorm(upper,mean=mean,sd=sd)),mean=mean,sd=sd)
}

n = 1000
a = 61965
D = KarnatakaForest$D
WD = myrtruncnorm(a, lower = 0.1, upper = 1.39, mean = 1)
coord12 = cbind( KarnatakaForest$long[1:a], KarnatakaForest$lat[1:a] )
D_simu = replicate(n, myrtruncnorm(a, mean= D[1:a], lower = 0.1, upper = 500))
WD_simu = replicate(n, myrtruncnorm(a, mean = WD[1:a], lower = 0.08, upper = 1.39))
selec <- sample(1:1001, n)






########## function coord
coord = function(n, coord, WD_simu, D_simu, selec){
  param_7 <- NULL
  data(param_7, envir = environment()) # posterior parameters from MCMC algorithm
  
  bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates
  
  # Posterior model parameters 
  intercept <- param_7[selec, "intercept"] # vector of simulated intercept
  coeffWD <- param_7[selec, "logwsg"] # vector of simulated parameter associated to ln(WD)
  coefflnD <- param_7[selec, "logdbh"] # vector of simulated parameter associated to ln(D) 
  coefflnD2 <- param_7[selec, "logdbh2"] # vector of simulated parameter associated to ln(D)^2
  coeffE <- -param_7[selec, "E"] # vector of simulated parameter associated to E
  coeffTmp <- param_7[selec, "temp"] # vector of of simulated parameter associated to tempsea coeff
  coeffCWD <- param_7[selec, "cwd"] # vector of of simulated parameter associated to CWD coeff
  coeffPS <- param_7[selec, "prec"] # vector of of simulated parameter associated to precSeas coeff
  RSE <- param_7[selec,"sd"] # vector of simulated RSE values
  
  # Recalculating n E values based on posterior parameters associated with the bioclimatic variables
  Tmp <- replicate(n, bioclimParams$tempSeas)
  CWD <- replicate(n, bioclimParams$CWD)
  PS <- replicate(n, bioclimParams$precSeas)
  
  Esim <- sweep(Tmp, MARGIN = 2, coeffTmp, "*") + sweep(CWD, MARGIN = 2, coeffCWD, "*") +
    sweep(PS, MARGIN = 2, coeffPS, "*")
  
  # Applying AGB formula over simulated matrices and vectors
  AGB_simu <- sweep(sweep(log(WD_simu), MARGIN = 2, coeffWD, "*") + 
                      sweep(log(D_simu), MARGIN = 2, coefflnD, "*") + 
                      sweep(log(D_simu)^2, MARGIN = 2, coefflnD2, "*")+
                      sweep(Esim, MARGIN = 2, coeffE, "*"),
                    MARGIN = 2, intercept, '+')
  return(AGB_simu)
}

# coord1 = function(n, coord, WD_simu, D_simu, selec){
#   param_7 <- NULL
#   data(param_7, envir = environment()) # posterior parameters from MCMC algorithm
#   
#   bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates
#   
#   AGB_simu = coord_cpp(tempSeas = bioclimParams$tempSeas, coeffTmp = param_7$temp[selec], 
#                        CWD = bioclimParams$CWD, coeffCWD = param_7$cwd[selec],
#                        PS = bioclimParams$precSeas, coeffPS = param_7$prec[selec],
#                        WD_simu = WD_simu, D_simu = D_simu, coeffWD =  param_7$logwsg[selec], coefflnD = param_7$logdbh[selec],
#                        coefflnD2 = param_7$logdbh2[selec], coeffE = param_7$E[selec], intercept = param_7$intercept[selec])
#   
#   return(AGB_simu)
#   
# }

coord1 = function(n, coord, WD_simu, D_simu, selec){
  param_7 <- NULL
  data(param_7, envir = environment()) # posterior parameters from MCMC algorithm
  
  bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates
  
  # Posterior model parameters 
  RSE <- param_7[selec,"sd"] # vector of simulated RSE values
  
  # Recalculating n E values based on posterior parameters associated with the bioclimatic variables
  Tmp <- replicate(n, bioclimParams$tempSeas)
  CWD <- replicate(n, bioclimParams$CWD)
  PS <- replicate(n, bioclimParams$precSeas)
  
  Esim <- t(Tmp) * param_7[selec, "temp"] + t(CWD) * param_7[selec, "cwd"] + t(PS) * param_7[selec, "prec"]
  
  # Applying AGB formula over simulated matrices and vectors
  logD = t(log(D_simu))
  AGB_simu <- t( t(log(WD_simu)) * param_7[selec, "logwsg"] + 
                   t(log(D_simu)) * param_7[selec, "logdbh"] + 
                   t(log(D_simu)^2) * param_7[selec, "logdbh2"] + 
                   Esim * -param_7[selec, "E"] + 
                   param_7[selec, "intercept"] )
  return(AGB_simu)
}


coord4 = function(n, coord, WD_simu, D_simu, selec){
  param_7 <- NULL
  data(param_7, envir = environment()) # posterior parameters from MCMC algorithm
  
  bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates
  
  # Posterior model parameters 
  RSE <- param_7[selec,"sd"] # vector of simulated RSE values
  
  # Recalculating n E values based on posterior parameters associated with the bioclimatic variables
  Esim <- t(bioclimParams$tempSeas %o% param_7[selec, "temp"] + 
    bioclimParams$CWD %o% param_7[selec, "cwd"] + 
    bioclimParams$precSeas %o% param_7[selec, "prec"])
  
  # Applying AGB formula over simulated matrices and vectors
  logD = t(log(D_simu))
  AGB_simu <- t( t(log(WD_simu)) * param_7[selec, "logwsg"] + 
                   t(log(D_simu)) * param_7[selec, "logdbh"] + 
                   t(log(D_simu)^2) * param_7[selec, "logdbh2"] + 
                   Esim * -param_7[selec, "E"] + 
                   param_7[selec, "intercept"] )
  return(AGB_simu)
}










#### Try of a new function coord

coord2 = function(n, coord, WD_simu, D_simu, selec){
  param_7 <- NULL
  data(param_7, envir = environment()) # posterior parameters from MCMC algorithm
  rownames(param_7) = NULL
  
  bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates
  
  # Posterior model parameters 
  RSE <- param_7[selec,"sd"] # vector of simulated RSE values
  
  # Recalculating n E values based on posterior parameters associated with the bioclimatic variables
  Esim = as.matrix( param_7[selec, c("temp", "prec", "cwd")] ) %*% t(as.matrix(bioclimParams))
  rownames(Esim) = NULL
  
  # Applying AGB formula over simulated matrices and vectors
  AGB_simu <- t( t(log(WD_simu)) * param_7[selec, "logwsg"] + 
                   t(log(D_simu)) * param_7[selec, "logdbh"] + 
                   t(log(D_simu)^2) * param_7[selec, "logdbh2"] + 
                   Esim * -param_7[selec, "E"] + 
                   param_7[selec, "intercept"] )
  return(AGB_simu)
  
}

coord3 = function(n, coord, WD_simu, D_simu, selec){
  param_7 <- NULL
  data(param_7, envir = environment()) # posterior parameters from MCMC algorithm
  rownames(param_7) = NULL
  
  bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates
  
  # Posterior model parameters 
  RSE <- param_7[selec,"sd"] # vector of simulated RSE values
  
  # Recalculating n E values based on posterior parameters associated with the bioclimatic variables
  Esim = tcrossprod(as.matrix(param_7[selec, c("temp", "prec", "cwd")]), as.matrix(bioclimParams))
  
  # Applying AGB formula over simulated matrices and vectors
  AGB_simu <- t( t(log(WD_simu)) * param_7[selec, "logwsg"] + 
                   t(log(D_simu)) * param_7[selec, "logdbh"] + 
                   t(log(D_simu)^2) * param_7[selec, "logdbh2"] + 
                   Esim * -param_7[selec, "E"] + 
                   param_7[selec, "intercept"] )
  return(AGB_simu)
  
}


all.equal( coord3(n, coord12, WD_simu, D_simu, selec), coord(n, coord12, WD_simu, D_simu, selec), 
           tolerance = 10^-14,
           check.attributes = F)


res = microbenchmark("original"= coord(n, coord12, WD_simu, D_simu, selec),
                     "modified1" = coord1(n, coord12, WD_simu, D_simu, selec),
                     "modified2" = coord2(n, coord12, WD_simu, D_simu, selec),
                     "modified3" = coord3(n, coord12, WD_simu, D_simu, selec),
                     "modified4" = coord4(n, coord12, WD_simu, D_simu, selec), times = 10)
