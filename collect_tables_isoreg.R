#############################################
##  Bootstrapping Isotonic regression Estimator
##  Cattaneo Jansson Nagasawa
##  mar-28-2023
##  Constructing tables
#############################################
rm(list=ls())
library(Hmisc)

models <- 1:3
n <- 1000

theta0 <- c(2,0,24)

n.outrow <- 7
out <- matrix(NA, nrow = n.outrow, ncol = 4*length(models))

for (m in models){
    results <- read.table( paste0("output/output_isoreg_m",m,".txt") )
    S <- nrow(results)/n.outrow
    subs <- c( ceiling(n^(1/2) ), ceiling(n^(2/3) ), ceiling(n^(4/5) )  )
    if (m==1|m==2) rate.adj <- rep(c(1,(subs/n)^(1/3),rep(1,n.outrow-4) ),S)
    if (m==3) rate.adj <- rep(c(1,(subs/n)^(3/7),rep(1,n.outrow-4) ),S)
    
    ci.low <- results[,"thetahat"] - results[,"q.975"]*rate.adj
    ci.upp <- results[,"thetahat"] - results[,"q.025"]*rate.adj
    test.tmp <- ( ci.low <= theta0[m]) & (theta0[m] <= ci.upp )
    ci.tmp <- ci.upp - ci.low
    
    ## average tuning parameters
    out[,4*m-3] <- by(results[,"D1"], results[,"label"], mean)
    out[,4*m-2] <- by(results[,"D3"], results[,"label"], mean)
    out[,4*m-1] <- by(test.tmp, results[,"label"], mean)
    out[,4*m] <- by(ci.tmp, results[,"label"], mean)
}    
    
output.N1 <- latex(round(out,3), file = "table_isoreg.txt",
                   append=FALSE, table.env=FALSE, center="none", title="",
                   n.cgroup=rep(4,length(models)),
                   cgroup=paste0("DGP ",models),
                   colheads= rep(c("$\\tilde{\\mathcal{D}}_{1,n}$","$\\tilde{\\mathcal{D}}_{3,n}$","Coverage", "Length"),length(models)),
                   n.rgroup = c(1,3,n.outrow-4),
                   rgroup = c("Standard","m-out-of-n","Reshaped"),
                   rowname=c("", "$m = \\lceil n^{1/2} \\rceil$", "$m = \\lceil n^{2/3} \\rceil$",
                             "$m = \\lceil n^{4/5} \\rceil$",
                             "Oracle", "ND known q", "ND robust")
)

