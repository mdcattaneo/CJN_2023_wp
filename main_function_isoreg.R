#############################################
##  Bootstrapping Isotonic Regression Estimator
##  Cattaneo Jansson Nagasawa
##  Mar-28-2023
##  Function File
#############################################
library(fdrtool)
library(lpdensity)
library(Rcpp)
sourceCpp("isoreg.cpp")
## setting seed
SeedFunction <- function() {
    args = commandArgs(trailingOnly=TRUE)
    
    # error handling
    if (length(args)==0) { args <- list("1") }
    
    # the argument is
    cat(paste("The argument is", args[[1]], "\n", sep=" "))
    
    # generate seed sequence
    set.seed(42)
    Seeds <- sample(1:50000, 50000, replace=FALSE)
    
    # set the seed
    cat(paste("The current seed is", toString(Seeds[as.integer(args[[1]])]), "\n", sep=" "))
    set.seed(Seeds[as.integer(args[[1]])])
    
}



##Data Generating Process
dgp <- function(N,model){
  A <- runif(N)
  if(model==1){
    Y <- 2*exp(A-0.5)  + rnorm(N)
  }
  if (model==2){
    Y <- 2*(A-0.5) + exp(A)*rnorm(N)
    
  }
  if (model==3){
    Y <- 24*( exp(A-0.5) - (A-0.5) - 0.5*(A-0.5)^2 ) + rnorm(N,sd=0.1)
  }
  out <- data.frame(y=Y,a=A)
  return(out)
}


## true constants of M0
D0 <- matrix(c(1,0,
               1,0,
               0,1), nrow=2,ncol=3)


## Compute isotonic regression estimate
Iso.reg <- function(A,PSI,x0){
    tmp <- isoreg(A,PSI)
    grd <- A[tmp$ord]
    ind <- which.max(grd > x0s[m] )-1
    if (ind>0){
      out <- tmp$yf[ind]
    }else{
      out <- out <- tail(tmp$yf,n=1)
    } 
    return(  out )
}


## Rule-of-thumb step size
ROT.step <- function(){
  obs$a0 <- obs$a -x0s[m]
  ref.lm <- lm(y~a0 +I(a0^2)+I(a0^3)+I(a0^4)+I(a0^5),data = obs)
  ref.gamma <- ref.lm$coefficients
  muAhat <- mean(obs$a)
  sigmaAhat <- sd(obs$a)
  dmu <- ref.gamma[2]
  d2mu <- 2*ref.gamma[3] 
  d3mu <- 6*ref.gamma[4] 
  d4mu <-24*ref.gamma[5]
  d5mu <- 120*ref.gamma[6]
  
  f <- dnorm((x0s[m]-muAhat)/sigmaAhat)/sigmaAhat
  df <- -(x0s[m]-muAhat)/sigmaAhat^3*dnorm((x0s[m]-muAhat)/sigmaAhat)
  d2f <- ( (x0s[m]-muAhat)^2/sigmaAhat^5 -1/sigmaAhat^3 )*dnorm((x0s[m]-muAhat)/sigmaAhat)
  d3f <- ( 3*(x0s[m]-muAhat)/sigmaAhat^5 -(x0s[m]-muAhat)^3/sigmaAhat^7 )*dnorm((x0s[m]-muAhat)/sigmaAhat)
  d4f <- ( 3/sigmaAhat^5 - 6*(x0s[m]-muAhat)^2/sigmaAhat^7 + (x0s[m]-muAhat)^4/sigmaAhat^9 )*dnorm((x0s[m]-muAhat)/sigmaAhat)
  
  ## j=1
  B1 <- -4/720*( 5*dmu*d4f + 10*d2mu*d3f+10*d3mu*d2f+5*d4mu*df+d5mu*f )
  V1 <- 113/144*(summary(ref.lm)$sigma)^2*f
  eps1 <- (3/8*V1/B1^2)^(1/11)*n^(-1/11)
  
  ## j=3
  B3 <- 5/720*(10*d3mu*d2f + 5*d4mu*df + d5mu*f)
  V3 <- 7/144*(summary(ref.lm)$sigma)^2*f
  eps3 <- (7/4*V3/B3^2)^(1/11)*n^(-1/11)
  
  return(c(eps1,eps3))
}



## numerical derivative estimator for M
ND.mhat1 <- function(eps){
  lambdas1 <- c(2/3,2/3,-1/24,-1/24)
  Y0 <- upsilon_hat(obs$a, obs$y, x0s[m],thetahat)
  Y1 <- upsilon_hat(obs$a, obs$y, x0s[m]+1*eps,thetahat) - Y0
  Y2 <- upsilon_hat(obs$a, obs$y, x0s[m]+(-1)*eps,thetahat) - Y0
  Y3 <- upsilon_hat(obs$a, obs$y, x0s[m]+2*eps,thetahat) - Y0
  Y4 <- upsilon_hat(obs$a, obs$y, x0s[m]+(-2)*eps,thetahat) - Y0
  max(c( eps^(-2)*( lambdas1%*%c(Y1,Y2,Y3,Y4) ) ),0)
}

ND.mhat3 <- function(eps){
  lambdas3 <- c(-1/6,-1/6,1/24,1/24)
  Y0 <- upsilon_hat(obs$a, obs$y, x0s[m],thetahat)
  Y1 <- upsilon_hat(obs$a, obs$y, x0s[m]+1*eps,thetahat) - Y0
  Y2 <- upsilon_hat(obs$a, obs$y, x0s[m]+(-1)*eps,thetahat) - Y0
  Y3 <- upsilon_hat(obs$a, obs$y, x0s[m]+2*eps,thetahat) - Y0
  Y4 <- upsilon_hat(obs$a, obs$y, x0s[m]+(-2)*eps,thetahat) - Y0
  max(c( eps^(-4)*( lambdas3%*%c(Y1,Y2,Y3,Y4) ) ),0)
}


#########################################################################################
## MERGE OUTPUT FILES
#########################################################################################
merge.output = function(cpus=100,models=1:3){
  
    for (m in models) {
        files = NULL;
        for (s in 1:cpus) {
            filename = paste0("output/parts/output_isoreg_m",m,"_",s,".txt")
            if (file.exists(filename)) files = c(files, filename)
        }
        if ( length(filename) > 0 ){
            write.table(do.call("rbind", lapply(files, function(x) read.table(x))), file = paste0("output/output_isoreg_m",m,".txt") )
            message("Model ",m," done -- ",length(files)," Files read.")
        } else message("Model ",m," not available.")
    }
}
