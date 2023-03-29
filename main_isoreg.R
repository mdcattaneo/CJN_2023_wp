#############################################
##  Bootstrapping Isotonic Regression Estimator
##  Cattaneo Jansson Nagasawa
##  Mar-28-2023
##  Main File
#############################################
rm(list=ls())

#########################################################################################
## MONTECARLO SETUP
#########################################################################################
ncpus <- 100

## sample size
n <- 1000
## Number of bootstrap replications
B <- 2000
## Total number of Monet Carlo iterations
S <- 4000%/%ncpus
## evaluation point of regression function: f(x0)
x0s <- rep(0.5,3)

ags <- commandArgs(trailingOnly = TRUE)
if (length(ags)<=1)ags <- c(1,1)
ag <- ags[1]

models <- as.integer(ags[2])

source("main_function_isoreg.R")
## SEED SETUP FOR PARALLELIZATION
SeedFunction()

num.rows.table <- 7 # 1 (std nonpar boot) + 3 (m-out-of-n) + 1 (Oracle) +  2 (ND: known q + robust)


col.names <- c("q.025","q.975","D1","D3","thetahat","label") 
#############################
## Monte Carlo Simulations ##
#############################
time <- Sys.time()

for (m in models){
    cat("\n##########################################")
    cat("\n###         Starting Model ",m,"         ###")
    cat("\n##########################################")
    
    ## matrix to store results
    out <- matrix(NA, nrow = S*num.rows.table, ncol = 6, dimnames = list(NULL, col.names) )
    eps.out <- rep(NA,S)
    Row <- 1
    ## subsample size for m-out-of-n bootstrap
    subsample <- c( ceiling(n^(1/2) ), ceiling(n^(2/3) ), ceiling(n^(4/5) )  )
    Morac <- D0[,m]
    
    for (s in 1:S){
        obs <- dgp(n,m)
        ## isotonic density estimate from the original sample
        thetahat <- Iso.reg(obs$a,obs$y,x0s[m])
        
        ## Numerical derivative estimator for the mean function
        # Rule-of-thumb step size
        ROT.eps <- ROT.step()
        
        M1hat <- ND.mhat1(ROT.eps[1])
        M3hat <- ND.mhat3(ROT.eps[2])
        
        # mean function estimate with known q
        if (m==1|m==2){
          Mhat.knowq <- c(M1hat,0)
        }
        if (m==3){
          Mhat.knowq <- c(0,M3hat)
        }
        
        # "robust" inference
        Mhat.robust <- c(M1hat,M3hat)
        
        ##########################
        ## Bootstrap Iterations ##
        ##########################
        # matrices to store the bootstrap computations
        std.boot <- rep(NA, B)
        mn.boot <- matrix(NA, nrow = B, ncol = length(subsample))
        oracle <- rep(NA, B)
        ND.knowq <- rep(NA, B)
        ND.robust <- rep(NA, B)
        
        for (b in 1:B){
            boot.obs <- obs[sample(n, replace = T),]
            
            ## standard nonparametric bootstrap
            std.boot[b] <- Iso.reg(boot.obs$a,boot.obs$y,x0s[m])
            
            ## m-out-of-n bootstrap
            for (mn in 1:length(subsample) ){
                mn.boot[b,mn] <- Iso.reg(boot.obs[1:subsample[mn],"a"],boot.obs[1:subsample[mn],"y"],x0s[m])
            }
            ## oracle mean function
            oracle[b] <- gridsearch(obs$a, obs$y, boot.obs$a,boot.obs$y,Morac,x0s[m],thetahat)
            
            ## ND
            ND.knowq[b] <- gridsearch(obs$a, obs$y, boot.obs$a,boot.obs$y,Mhat.knowq,x0s[m],thetahat)
            ND.robust[b] <- gridsearch(obs$a, obs$y, boot.obs$a,boot.obs$y,Mhat.robust,x0s[m],thetahat)
            
        }
        ## end of bootstrap iterations
        
        ##compute quantiles
        for (k in 1:num.rows.table){
            if (k == 1){ #standard nonpara bootstrap
                out[Row, 1:2] <- quantile(std.boot - thetahat, probs = c(0.025,0.975), type = 1)
                out[Row, 5] <- thetahat
                out[Row, 6] <- k
            } 
            if (2 <= k & k <= 4){ #m-out-of-n bootstrap
                out[Row, 1:2] <- quantile( mn.boot[,k-1] - thetahat, probs = c(0.025, 0.975), type = 1 )
                out[Row, 5] <- thetahat
                out[Row, 6] <- k
            }
            if (k==5){ # Bootstrap with oracle M0
               out[Row,1:2] <- quantile( oracle - thetahat, probs = c(0.025, 0.975), type = 1 )
               out[Row, 3:4] <- Morac
               out[Row, 5] <- thetahat
               out[Row, 6] <- k
            }
            if (k==6){ #ND known q
              out[Row, 1:2] <- quantile(ND.knowq-thetahat, probs = c(0.025, 0.975), type = 1 )
              out[Row, 3:4] <- Mhat.knowq
              out[Row, 5] <- thetahat
              out[Row, 6] <- k
            }
            if (k==7){ #ND robust
              out[Row, 1:2] <- quantile(ND.robust-thetahat, probs = c(0.025, 0.975), type = 1 )
              out[Row, 3:4] <- Mhat.robust
              out[Row, 5] <- thetahat
              out[Row, 6] <- k
            }
            Row = Row + 1
        }
        cat(paste("\nSimulations Completed:",s,
                  "Diff:", round(difftime(Sys.time(),time,units="mins"),2),"mins"));
        time <- Sys.time();
    }
    
    ##################
    ## Save Results ##
    ##################
    write.table(out, file = paste0("output/parts/output_isoreg_m",m,"_",ag,".txt") )
}

merge.output(ncpus, models)
