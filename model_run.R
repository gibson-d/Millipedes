######################################################################
# Section 2: Model likelihoods 
######################################################################
library(nimble)
# Specify model in NIMBLE language
model_1 <- nimbleCode({
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Specify priors
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      mean.gamma ~ dbeta(1,1)
        mean.phi ~ dbeta(1,1)
        
      mu.phi <- logit(mean.phi)
      mu.gamma <- logit(mean.gamma)
 
      beta.temp[1] ~ dnorm(0, 1)
      beta.temp[2] ~ dnorm(0, 1)
      
       for (t in 1:(n.occasions-1)) {
             logit(phi[t])  <- mu.phi + beta.temp[1] * temp[t]
          logit(gammap[t])  <- mu.gamma + beta.temp[2] * temp[t]
     }
        psi  ~ dbeta(1,1)
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
     # Secondary occasions p's
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        mean.p ~ dbeta(1,1)
         sig.p ~ dunif(0,5)
        
    shape <- pow(sig.p,-2) * mean.p
     rate <- pow(sig.p,-2) * (1 - mean.p)
     
    for (t in 1:n.occasions){
            lp[t] ~ dbeta(shape,rate)
      for (k in 1:n.sec[t]){
           p[t,k] <- lp[t] 
        hits[t,k] ~ dbin(p[t,k], hits_misses[t,k])
        }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      # Primary occasions p's  
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
       pstar[t] <- 1 - ( (1 - p[t,1]) * (1 - p[t,2]) )
      }
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
     # Prior for entry probabilities
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        for (t in 1:n.occasions) {
              beta[t]  ~ dgamma(1, 1)
                 b[t] <- beta[t] / sum(beta[1:n.occasions])
        }
         eta[1] <- b[1] 
        for (t in 2:n.occasions) {
         eta[t] <- b[t] / (1 - sum(b[1:(t-1)]))
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ########### State Process ############
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        for (i in 1:M){                                                   # Loop through real + psuedo individuals
         w[i] ~ dbern(psi)                                                # Latent inclusion state
         z[i,1] ~ dbern(eta[1])                                           # Entry into first occasion
        
        for (t in 2:n.occasions) {
          q[i,t-1] <- 1 - z[i,t-1]                                        # Have not entered yet 
          mu1[i,t] <- phi[t-1] * z[i,t-1] + eta[t] * prod(q[i,1:(t-1)]) # Presence = survived since last + entered since last
          z[i,t] ~ dbern(mu1[i,t])
          }
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ######### Observation Process ########
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        for (i in 1:n.ch){                                            # Loop through real individuals
          y[i,1] ~ dbern(z[i,1] * pstar[1] * w[i])                    # Latent state * detection * latent inclusion 
        for (t in 2:n.occasions){
         y[i,t] ~ dbern(z[i,t] * pstar[t] * w[i] * gammap[t-1])      # Latent state * detection * latent inclusion * availability 
        }
       }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ########## Chandler trick ############
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        for(i in (n.ch+1):M) {                                      # Loop through psuedo individuals
          pzo[i,1] <- 1 - (pstar[1] * z[i,1] * w[i])                # probability of not being detected on first occasion
        for (t in 2:n.occasions){
          pzo[i,t] <- 1 - (pstar[t] * z[i,t] * w[i] * gammap[t-1])  # probability of not being detected on each other occasion
        }
          zeros[i] ~ dbern(1-prod(pzo[i,1:n.occasions]))       # probability of at being detected at least once (zeros = vector of 0's for each augmented individual)
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        ################# Calculate derived population parameters ##################
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
          for (i in 1:M){
             for (t in 1:n.occasions){
                u[i,t] <- z[i,t]*w[i]     # Deflate latent state (u)
              zstop[i] <- sum(u[i,1:(n.occasions-1)])    # individual stopover duration (known + pseudo)
                }
             }
               zes <- mean(zstop[1:nch]                  # mean stopover, known-only
          for (i in 1:M){
             recruit[i,1] <- u[i,1]
             for (t in 2:n.occasions){
                recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
                } #t
             } #i
          for (t in 1:n.occasions){
             N[t] <- sum(u[1:M,t])        # Actual population size
             B[t] <- sum(recruit[1:M,t])  # Number of entries
             } #t
          for (i in 1:M){
               Nind[i] <- sum(u[i,1:n.occasions])
             Nalive[i] <- 1-equals(Nind[i], 0)
             } #i
          Nsuper <- sum(Nalive[1:M])         # Superpopulation size
    

  
    }#i
    for (i in (loop.limits[j] + 1): (max(loop.limits)+20)){
    zstop[i,j] <- 99999                                            # 'ragged edged objects'protection
    }#i
    zes[j] <- mean(zstop[1:loop.limits[j],j])
    }    #j
    
    
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      ############################## End model ################################### 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
          
    })

ch[ch == 2] <- 0 
ch.js <- ch.ms                             #Transistion MS-ch to JS-ch
ch.js[ch.js == 2] <- 0                     #Augmented capture history (JS_S format)
ch.js <- ch.js[,-1]                        
######################################################################
z.init <-ch.js                             # Set initial values for latent state     
z.init[z.init==0] <- 1
w.init <- rep(1, nrow(ch.js))              # Set initial values for latent inclusion state

######################################################################
# Bundle data
dat       <- list(y = ch,  hits = alpha, hits_misses = (alpha + beta), zeros = c(rep(NA, nrow(ch)), rep(0, nz)),lengths = c(size, rep(NA, nz)), temp  = temp.min[,1] )
constants <- list( n.occasions = dim(ch.js)[2], M = dim(ch.js)[1],   n.sec = rep(2,13), n.ch = dim(ch)[1] )

pars <- c("b", "beta", "mean.phi", "mean.gamma", "psi", 'p', 'pstar','Nsuper','N', 'beta.temp','zes', 'zstop')
# Initial values
 initsFunction <- function()list(z = z.init, w = w.init, mean.gamma = 0.3, 
                                 mean.phi = 0.8, beta.temp = c(0,0),
                                 mean.p = 0.25, sig.p = .25, 
                                 psi = 0.5) 

inits <- initsFunction() 
start_time <- Sys.time()

library(parallel)
library(coda)
nc      <- 3 # number of chains

cl_new5 <-makeCluster(nc,timeout=5184000)
clusterExport(cl_new5, c("model_1", "inits", "dat", "constants", "pars"))

for (j in seq_along(cl_new5)) {
  set.seed(j)
  inits <- initsFunction() 
  clusterExport(cl_new5[j], "inits")
}
out <- clusterEvalQ(cl_new5, {
  library(nimble)
  library(coda)
  
  model <- nimbleModel( code = model_1, constants = constants,  dat =  dat, inits = inits)
  
  model$simulate(c('p','lp', 'pstar', 'beta', 'b', 'eta', 'mu1','pzo','zeros'))
  model$initializeInfo()
  model$calculate()
  
  mcmc_Conf  <- configureMCMC(model, useConjugacy = FALSE, thin2 = 5, 
                              monitors2 = c('pstar','N','B','Nsuper','mean.phi','mean.gamma', 'phi', 'gammap','b','beta', 'shape','rate'))
  modelMCMC  <- buildMCMC( mcmc_Conf)
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC)
  
  CmodelMCMC$run(100000, thin = 1, thin2 = 1, nburnin = 50000)
  
  return(list( as.mcmc(as.matrix(CmodelMCMC$mvSamples )),
               as.mcmc(as.matrix(CmodelMCMC$mvSamples2))))
})

end_time <- Sys.time()
end_time - start_time

samples2 <- list( chain1 =   out[[1]][[1]], chain2 =   out[[2]][[1]], chain3 =   out[[3]][[1]])
samples1 <- list( chain1 =   out[[1]][[2]], chain2 =   out[[2]][[2]], chain3 =   out[[3]][[2]])

mcmcList1 <- as.mcmc.list(lapply(samples1, mcmc))
mcmcList2 <- as.mcmc.list(lapply(samples2, mcmc))
