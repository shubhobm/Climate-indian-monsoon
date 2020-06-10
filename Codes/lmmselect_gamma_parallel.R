## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
#setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/Rain")
# setwd('D:/Study/My projects/Climate-indian-monsoon/Codes')
source('misc_functions.R')

library(rrcov)
library(fda.usc)
library(lme4)
# library(MuMIn)
library(parallel)

# read in data
# rainsmall = read.csv("../data/rainsmall.csv")
rainsmall = read.csv("../Data/rainsmall.csv")

# check full model
rainsmall[-(1:3)] = scale(rainsmall[-(1:3)])
trainset = which(rainsmall$year<2003)
testset = which(rainsmall$year>=2003)
varnames = names(rainsmall)[-(1:3)]
formula = paste(varnames, collapse="+")
random_terms = "+ (1|year)"
formula = as.formula(paste("log(PRCP+1) ~", formula, random_terms))
mod.full = lmer(formula, data=rainsmall, subset=trainset)
# summary(mod.full)
# anova(mod.full)
# r.squaredGLMM(mod.full)

y = getME(mod.full, 'y')
x = getME(mod.full, 'X')
n = nrow(x)
p = ncol(x)-1

get.err = function(sdn){
  require(lme4)
  require(fda.usc)
  nboot = 1e3
  set.seed(07152015)
  sdn = sdn

  ## loop to get full and drop-1 bootstrap estimates
  loopfun = function(b){
    require(lme4)
    
    # get model corresponding to bootstrap replicate
    set.seed(b)
    assign("w", rgamma(length(trainset),sdn,sdn), envir=.GlobalEnv)
    bmod = lmer(formula, data=rainsmall[trainset,], weights=w)
    
    # store model estimates
    beta.mat = matrix(0, nrow=p+1, ncol=p)
    beta.hat.b = getME(bmod, 'beta')[-1]
    beta.mat[p+1,] = beta.hat.b
    for(j in 1:p){
      beta.mat[j,-j] = beta.hat.b[-j]
    }
    beta.mat
  }
  
  # run function
  output = lapply(1:nboot, loopfun)
  
  # process output into several matrices, corresponding to full and drop-1 model estimates
  beta.mat.full = matrix(0, nrow=nboot, ncol=p)
  jbeta.list = list()
  for(j in 1:p){
    jbeta.list[[j]] = beta.mat.full
  }
  
  for(b in 1:nboot){
    output.b = output[[b]]
    beta.mat.full[b,] = output.b[p+1,]
    for(j in 1:p){
      jbeta.list[[j]][b,] = output.b[j,]
    }
  }
  
  # now calculate depth
  SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
  SSPmat.d[,p+1] = mdepth.RP(beta.mat.full, beta.mat.full)$dep
  for(j in 1:p){
    SSPmat.d[,j] = mdepth.RP(jbeta.list[[j]], beta.mat.full)$dep
  }
  
  # get p-values
  pVal = rep(1, p+1)
  for(i in 1:p){
    pVal[i] = t.test(SSPmat.d[,i], SSPmat.d[,i+1], paired=TRUE)$p.value
  }
  
  Cn.frame = data.frame(DroppedVar = c(paste("-", names(data.frame(x))[-1]), "<none>"),
                        Cn = apply(SSPmat.d, 2, mean))
  Cn.frame = Cn.frame[with(Cn.frame, order(Cn)),]
  row.names(Cn.frame) = NULL
  
  # final model, and return prediction error
  noneCn = Cn.frame$Cn[which(Cn.frame$DroppedVar == "<none>")]
  which.final = which(apply(SSPmat.d, 2, mean) < .95*noneCn)
  # Cn.frame$ind = "N"
  # Cn.frame$ind[which(Cn.frame$ind < .95*noneCn)] = "Y"
  fixed.final = paste(varnames[which.final], collapse="+")
  form.final = as.formula(paste("log(PRCP+1) ~", fixed.final, random_terms))
  mod.final = lmer(form.final, data=rainsmall, subset=testset)
  list(Cn.frame,
       mod.final,
       mean((y - x[,which.final] %*% getME(mod.final,'beta')[-1])^2))
}

# get errors for a range of sdn values
sdn.vec = n^seq(.01,.16,by=.01)
outputs = mclapply(sdn.vec, get.err, mc.cores=16)
best.index = which.min(lapply(outputs, function(x) x[[3]]))
save(outputs[[best.index]], file="bestmodel.Rda")
outputs[[best.index]][[1]]
