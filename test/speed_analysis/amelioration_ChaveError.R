library(microbenchmark)

D_simu = replicate(1000, rnorm(50000, mean=15, sd=14.5))



myrtruncnorm <- function(n,lower = -1, upper = 1,mean=0,sd=1) {
  qnorm(runif(n,pnorm(lower,mean=mean,sd=sd),pnorm(upper,mean=mean,sd=sd)),mean=mean,sd=sd)
}





rt_msm = function(D_simu){
  fivePercent <- round( nrow(D_simu) * 5 / 100 )
  chaveError <- function(x)
  { 
    ## Assigning large errors on 5% of the trees
    largeErrSample <- sample(length(x), fivePercent)
    x[largeErrSample] <- myrtruncnorm(n = fivePercent, mean = x[largeErrSample], sd = 4.64, lower = 0.1, upper = 500)
    ## Assigning small errors on the remaining 95% trees
    D_sd <- 0.0062 * x[-largeErrSample] + 0.0904
    x[-largeErrSample] <- msm::rtnorm(n =length(D_sd), mean = x[-largeErrSample], sd = D_sd, lower = 0.1, upper = 500)
    return(x)
  }
  D_simu <- apply(D_simu, 2, chaveError)
  return(D_simu)
}

rt_myrtroncnorm = function(D_simu){
  fivePercent <- round( nrow(D_simu) * 5 / 100 )
  chaveError <- function(x)
  { 
    ## Assigning large errors on 5% of the trees
    largeErrSample <- sample(length(x), fivePercent)
    x[largeErrSample] <- myrtruncnorm(n = fivePercent, mean = x[largeErrSample], sd = 4.64, lower = 0.1, upper = 500)
    ## Assigning small errors on the remaining 95% trees
    D_sd <- 0.0062 * x[-largeErrSample] + 0.0904
    x[-largeErrSample] <- myrtruncnorm(n =length(D_sd), mean = x[-largeErrSample], sd = D_sd, lower = 0.1, upper = 500)
    return(x)
  }
  D_simu <- apply(D_simu, 2, chaveError)
  return(D_simu)
}


rt_improvement = function(D_simu){
  fivePercent <- round( nrow(D_simu) * 5 / 100 )
  chaveError <- function(x)
  { 
    ## Assigning large errors on 5% of the trees
    largeErrSample <- sample(length(x), fivePercent)
    
    D_sd = 0.0062 * x + 0.0904 # Assigning small errors on the remaining 95% trees
    D_sd[largeErrSample] = 4.64
    
    x <- myrtruncnorm(n =length(x), mean = x, sd = D_sd, lower = 0.1, upper = 500)
    return(x)
  }
  D_simu <- apply(D_simu, 2, chaveError)
  return(D_simu)
}

result = microbenchmark("rt_msm" = rt_msm(D_simu), 
                        "rt_myr" = rt_myrtroncnorm(D_simu), 
                        "rt_improvement"=rt_improvement(D_simu), times = 10)
plot(result)






tnorm = function(n, mean = 0, sd = 1, lower = -1, upper = 1){
  out = rnorm(n, mean, sd)
  
  a = which(out < lower || out > upper)
  
  while( length( a ) != 0){
    out[a] = rnorm( length(a), mean = mean[a], sd[a] )
    a = which(out < lower || out > upper)
  } 
  
  return(out)
}


myrtruncnorm <- function(n,lower = -1, upper = 1,mean=0,sd=1) {
  qnorm(runif(n,pnorm(lower,mean=mean,sd=sd),pnorm(upper,mean=mean,sd=sd)),mean=mean,sd=sd)
}


microbenchmark("rtnorm" = msm::rtnorm(10^5, lower = -1, upper = 1), 
                    "tnorm"= tnorm(10^5), 
                    "mytruncnorm"=myrtruncnorm(10^5, a=-1, b=1),
                    "rtruncnorm" = truncnorm::rtruncnorm(10^5, a=-1, b=1))





rt_improvement = function(D){
  
  D_simu = replicate(1000, D)
  
  fivePercent <- round( nrow(D_simu) * 5 / 100 )
  chaveError <- function(x)
  { 
    ## Assigning large errors on 5% of the trees
    largeErrSample <- sample(length(x), fivePercent)
    
    D_sd = 0.0062 * x + 0.0904 # Assigning small errors on the remaining 95% trees
    D_sd[largeErrSample] = 4.64
    
    x <- myrtruncnorm(n =length(x), mean = x, sd = D_sd, lower = 0.1, upper = 500)
    return(x)
  }
  D_simu <- apply(D_simu, 2, chaveError)
  return(D_simu)
}

rt_i = function(D){
  fivePercent <- round( length(D) * 5 / 100 )
  chaveError <- function(x)
  { 
    ## Assigning large errors on 5% of the trees
    largeErrSample <- sample(length(x), fivePercent)
    
    D_sd = 0.0062 * x + 0.0904 # Assigning small errors on the remaining 95% trees
    D_sd[largeErrSample] = 4.64
    
    x <- myrtruncnorm(n =length(x), mean = x, sd = D_sd, lower = 0.1, upper = 500)
    return(x)
  }
  D_simu = replicate(1000, chaveError(D))
  
  return(D_simu)
}

rt_i(D)

res = microbenchmark("rt_improvement"=rt_improvement(D), "rt_i"=rt_i(D))
