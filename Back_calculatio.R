rm(list = ls())
library(nimble)
library(coda) # For manipulation of MCMC results.
library(mgcv)
library(ggplot2)
library(ggpubr)
#library(tidyverse) # not on the cluster
library(ggmcmc)
#library(ggpubr)
library(data.table)
library(reshape)
library(dplyr)


path<-file.path("~","Nunavut_data.csv")
my_data = read.table(path,header = TRUE, sep = ",",stringsAsFactors = FALSE)
nobs = length(my_data$new_cases)

lmod_code = nimbleCode({
  for(i in 1:n){
    I[i] ~ dpois(lambda[i])
    # p[i] <- 1 - (lambda[i]/variance[i])
    # r[i] <- (lambda[i] * (1 - p[i])) / p[i]
    # I[i] ~ dnegbin(p[i], r[i])
    for (j in 1:i)
    {
      Z[i, j] <- f[(i+1)-j] * I[j]
    }
    
    # mu1[i] <- sum(Z[i, ])
    # # var1[i] <- mu1[i]+1
    # var1[i] ~ T(dnorm(mu1[i], varsigma), mu1[i],)
    # p1[i] <- 1 - (mu1[i]/var1[i])
    # r1[i] <- (mu1[i] * (1 - p1[i])) / p1[i]
    D[i] ~ dpois(sum(Z[i,]))
    
    sigma[i+1] ~ T(dnorm(1, 1), 0,)
    lambda[i+1] ~ T(dnorm(lambda[i], sigma[i]), 0,) 
    # variance[i+1] ~ T(dnorm(lambda[i+1], varsigma), lambda[i],)
    
  }
})



my_data = my_data$new_cases
nobs = length(my_data)


prob_vec = rep(0, nobs)
for(i in 1:nobs){
  prob_vec[i] = dlnorm(i, 1.434, 0.661)
}
plot(prob_vec)


Zic = matrix(data = rep(0, nobs*nobs), nrow = nobs, ncol = nobs)
inits1=list(lambda = rep(1, nobs+1), sigma=rep(1, nobs+1),I = my_data)
inits2=list(lambda = rep(2, nobs+1),  sigma=rep(1, nobs+1),I = my_data)
inits3=list(lambda = rep(0.5, nobs+1),  sigma=rep(1, nobs+1),I = my_data)
inits4=list(lambda = rep(1.5, nobs+1), sigma=rep(1, nobs+1),I = my_data)
inits5=list(lambda = rep(2.5, nobs+1), sigma=rep(1, nobs+1),I = my_data)

inits=list(chain1=inits1, chain2=inits2, chain3=inits3, chain4=inits4, chain5=inits5)
lmod_constants = list(n = nobs, varsigma = 1, f = prob_vec, Z=Zic)
lmod_data = list(D = my_data)

model <- nimbleModel(lmod_code, lmod_constants, lmod_data, inits)
compiled_model <- compileNimble(model, resetFunctions = TRUE)
mcmc_conf <- configureMCMC(compiled_model, monitors=c("I"), print=T)
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = model)

mcmc_samples = runMCMC(compiled_mcmc, inits=inits, nchains=5, nburnin=10000, niter=20000, thin=1,
                       samplesAsCodaMCMC = TRUE, summary=T, setSeed=c(1, 2, 3, 4, 5))

df = data.table(do.call('rbind', mcmc_samples$samples))
colMeans(df) 
sum(colMeans(df))

write.table(df,file = "Nunavut_Infection_estimates.csv",sep = ",", row.names = FALSE, col.names = FALSE)


p<-ggplot(df, aes(x=as.factor(get("I[30]")))) + 
  geom_histogram(color="black", fill="white", stat="count")
p




lstr = c("I[20]", "I[21]", "I[22]", "I[23]")
trace_plots = df[, ..lstr]
chainsPlot(trace_plots,file = "myplot.pdf",width = 7)

# trace_plots = trace_plots %>% mutate(itr=rep(seq(1, 20000), 3)) %>% mutate(chain=rep(c(1, 2, 3, 4, 5), 1, each=4000))

