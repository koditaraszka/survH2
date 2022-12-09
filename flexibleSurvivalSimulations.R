library(argparse)
library(MASS)
library(coxmeg)
library(dplyr)
library(stats)
library(survival)

#rm(list=ls())
#gc()

#parsing function: contains argparse from python
parsing = function(){
  parser <- ArgumentParser(description='Survival Simulation Framework')
  parser$add_argument('-N','--sampleSize', dest='N', type="integer", default = 500, required = F, help='number of individuals (default is 500)')
  parser$add_argument('-M','--numSNPs', dest='M', type="integer", default = 100, required = F, help='number of SNPs (default is 100)')
  parser$add_argument('-k','--causalSNPs', dest='k', type="integer", default = 100, required = F, help='number of causal SNPs (default is 100)')
  parser$add_argument('-s','--sims', dest='sims', type="integer", default = 100, required = F, help='number of simulations (default is 100)')
  parser$add_argument('-h2','--herit', dest='h2', type="double", default = 0.5, required = F, help='heritability/variance explained (default is 0.5)')
  parser$add_argument('-c','--censoring', dest='censor', type="double", default = 0, required = F, help='proportion of patients censored (default is 0/none)')
  parser$add_argument('-d','--dist', dest='dist', type="character", nargs = "*", default = c('exponential', 1), required = F, help='distribution of baseline hazard; 
                    options are [exponential lambda, weibull a lambda, lognorm mu sigma] (default is exponential 1)')
  parser$add_argument('-g','--geneticEffect', dest='effect', type="character", default = 'none', required = F, help='genetic effect none, SUM of fixed effects, or drawn from multivariate normal; 
                    options are [none, SUM, MVN]')
  parser$add_argument('--catCov', dest='catCov', type="character", default = 'none', nargs = "*", required = F, help='Categorial covariate default none or 4 parts: 
                    number of categories, mean effect size, std deviation, proportion with positive effect') 
  parser$add_argument('--catBeta', dest='catBeta', type="character", default = 'none', nargs = "*", required = F, help='Categorial covariate betas: either none or 1 beta for each categories -1')
  parser$add_argument('--binCov', dest='binCov', type="character", default = 'none', nargs = "*", required = F, help='Binary covariate default none or 5 parts: 
                    number of covariates, proportion of carriers, mean effect size, std deviation, proportion with positive effect')
  parser$add_argument('--binBeta', dest='binBeta', type="character", default = 'none', nargs = "*", required = F, help='Binary covariate betas: either none or 1 beta for each indicator')
  parser$add_argument('--tdCov', dest='tdCov', type="character", default = 'none', nargs = "*", required = F, help='time dependent covariates (default none) TODO')
  parser$print_help()
  return(parser)
  
}

#checkArgs function: sanity check input to parser
#Arguments: either args or commandLine arguments
checkArgs = function(arguments){
  if(arguments$k!=100) stop("TODO: Currently number of causal SNPs needs to be 100 b/c it's not actually used yet")
  if(!(arguments$dist[1] %in% c('exponential', 'weibull', 'lognorm'))) stop("Currently, distributions (-d, --dist) are exponential, weibull, or lognorm")
  if(!(arguments$effect) %in% c('SUM', 'MVN', 'none')) stop ('Currently, the genetic effect (-g, --geneticEffect) is either sum of fixed effects or drawn from multivariate normal')
  if(arguments$tdCov!='none') stop('Currently tdCov has not been created yet!')
}

#random function: set the random genetic effect for frailty model
#Arguments: N=#individuals, M=#SNPs, genEffect =[fixed or random] (either sum of fixed effects or from multivariate normal), h2=variance explained/tau
#Returns: List with $genetic = genetic effect on rate, $GRM = genetic relatedness matrix
random = function(N, M, genEffect, h2){
  if(genEffect=='none') return(list('genetic'=NULL,'GRM'=NULL))
  MAF = runif(M, min=0.01, max=0.99)
  G = NULL
  for(m in 1:M){
    g = rbinom(N, 2, prob = MAF[m])
    G = cbind(G, g)
  }
  X = scale(G)
  GRM = X %*% t(X) / M
  if(genEffect == 'SUM'){
    beta = rnorm(M,0,1)
    xb = sqrt(h2)*scale(X %*% beta)
    return(list('genetic'=xb,'GRM'=GRM))
  } else{
    xb=sqrt(h2)*scale(mvrnorm(1, mu = rep(0,N), Sigma = GRM))
    return(list('genetic'=xb,'GRM'=GRM))
  }
}

#fixed function: set the fixed effect covariates for frailty model
#Arguments: N=#individuals, catCov=paramters for categorical covariate or 'none', catBeta=set Betas for categorical variable if not random draws or none, 
#           binCov=parameters for binary covariate or 'none', binBeta=set Betas for binCov if not random draws or none, tdCov=paramters for time dependent covariate or 'none'
#Returns: List with $betas = effect size for each covariate, $covars = matrix of covariates
fixed = function(N, catCov, catBeta, binCov, binBeta, tdCov){
  covariates = NULL
  betas = NULL
  if(catCov[1]!='none'){
    if(length(catCov)==1 & catBeta[1]=='none') stop('--catCov has 1 argument (not none) but --catBeta is none. Either pass 4 arguments to --catCov or 1 argument and use --catBeta')
    if(length(catCov)==4 & catBeta[1]!='none') stop('--catCov has 4 arguments and --catBeta is not none. You cannot use both. Either pass 4 arguments to --catCov or 1 argument and use --catBeta')
    if(length(catCov)!=1 & length(catCov)!=4) stop('Either pass 4 arguments to --catCov or 1 argument and use --catBeta')
    categories = as.numeric(catCov[1])
    if(catBeta[1]!='none'){
      if(length(catBeta)!=(categories-1)) stop('--catBeta is not none but the number of arguments does not equal # of categories -1')
      betas = c(betas, as.numeric(catBeta))
    } else{
      meanEffect = as.numeric(catCov[2])
      sdEffect = as.numeric(catCov[3])
      posDir = as.numeric(catCov[4])
      negDir = 1-posDir
      effect = rnorm((categories-1), meanEffect, sdEffect)
      direction = sample(c(-1,1), (categories-1), replace = T, prob = c(posDir, negDir))
      betas = c(betas, effect*direction)
    }
    prop = runif(categories)
    prop = prop/sum(prop)
    while(min(prop)<0.05){
      prop = runif(categories)
      prop = prop/sum(prop)
    }
    catVar = sample(1:categories, N, replace = T, prob = prop)
    dummy = matrix(0, nrow=N, ncol=(categories-1))
    for(i in 1:(categories-1)){
      dummy[which(catVar==i),i] = 1
    }
    covariates = cbind(covariates, dummy)
  }
  if(binCov[1]!='none'){
    if(length(binCov)==1 & binBeta[1]=='none') stop('--catCov has 1 argument (not none) but --binBeta is none. Either pass 5 arguements to --binCov or 1 arguement and use --binBeta')
    if(length(binCov)==5 & binBeta[1]!='none') stop('--catCov has 5 argument and --binBeta is not none. You cannot use both. Either pass 5 arguements to --binCov or 1 arguement and use --binBeta')
    if(length(binCov)!=2 & length(binCov)!=5) stop('Either pass 5 arguments to --binCov or 2 argument and use --binBeta')
    count = as.numeric(binCov[1])
    prop = as.numeric(binCov[2])
    if(binBeta[1]!='none'){
      if(length(binBeta)!=count) stop('--binBeta is not none but the number of arguments does not equal # of binary variables')
      betas = c(betas, as.numeric(binBeta))
    } else{
      meanEffect = as.numeric(binCov[3])
      sdEffect = as.numeric(binCov[4])
      posDir = as.numeric(catCov[5])
      negDir = 1-posDir
      effect = rnorm(count, meanEffect, sdEffect)
      direction = sample(c(-1,1), categories, replace = T, prob = c(posDir, negDir))
      betas = c(betas, binBeta)
    }
    for(i in 1:count){
      bern = rbinom(N, 1, prop)
      covariates = cbind(covariates, bern)
    }
  }
  return(list("betas"=betas,"covars"=covariates))
}

#exponential function: generate data according to the exponential distribution
#Arguments: N=#individuals, fixed=fixed effects, random=genetic effect, censor=censoring rate, lambda=parameter for exponential
#Returns: List with $y as the time to event, $death = indicator if event occurs
exponential = function(N, fixed, random, censor, lambda){
  myrates = NULL
  if(is.null(fixed$betas) & is.null(random$genetic)) stop('No covariates or genetic effect arguments')
  if(is.null(fixed$betas)){
    myrates = random$genetic
  } else{
    myrates = fixed$covar %*% fixed$betas
    if(!is.null(random$genetic)){
      myrates = myrates + random$genetic
    }
  }
  U = runif(N, 0, 1)
  y_surv = (-log(1-U) / (as.numeric(lambda) * exp(myrates)))
  if(censor!=0){
    cen <- rexp(N, rate = censor)
    y_cen <- pmin(y_surv, cen)
    death = as.numeric(y_surv <= cen)
    return(list('y'=y_cen, 'death'=death))
  } else return(list('y'=y_surv, 'death'=rep(1,N)))
}

#weibull function: generate data according to the weibull distribution
#Arguments: N=#individuals, fixed=fixed effects, random=genetic effect, censor=censoring rate, a & lamba = paramters for weibull
#Returns: List with $y as the time to event, $death = indicator if event occurs
weibull = function(N, fixed, random, censor, gamma, lambda){
  myrates = NULL
  if(is.null(fixed$betas) & is.null(random$genetic)) stop('No covariates or genetic effect arguments')
  if(is.null(fixed$betas)){
    myrates = random$genetic
  } else{
    myrates = fixed$covar %*% fixed$betas
    if(!is.null(random$genetic)){
      myrates = myrates + random$genetic
    }
  }
  U = runif(N, 0, 1)
  y_surv = (-(1/as.numeric(lambda))*exp(myrates)*log(1-U))^(1/as.numeric(gamma))
  if(censor!=0){
    cen <- rexp(N, rate = censor)
    y_cen <- pmin(y_surv, cen)
    death = as.numeric(y_surv <= cen)
    return(list('y'=y_cen, 'death'=death))
  } else return(list('y'=y_surv, 'death'=rep(1,N)))
}

#lognormal function: generate data according to the lognormal distribution
#Arguments: N=#individuals, fixed=fixed effects, random=genetic effect, censor=censoring rate, a & lamba = paramters for weibull
#Returns: List with $y as the time to event, $death = indicator if event occurs
lognormal = function(N, fixed, random, censor, mu, sigma){
  myrates = NULL
  if(is.null(fixed$betas) & is.null(random$genetic)) stop('No covariates or genetic effect arguments')
  if(is.null(fixed$betas)){
    myrates = random$genetic
  } else{
    myrates = fixed$covar %*% fixed$betas
    if(!is.null(random$genetic)){
      myrates = myrates + random$genetic
    }
  }
  U = runif(N, 0, 1)
  y_surv = exp(as.numeric(mu) + as.numeric(sigma)*qnorm(1-exp((log(1-U)/exp(myrates)))))
  if(censor!=0){
    cen <- rexp(N, rate = censor)
    y_cen <- pmin(y_surv, cen)
    death = as.numeric(y_surv <= cen)
    return(list('y'=y_cen, 'death'=death))
  } else return(list('y'=y_surv, 'death'=rep(1,N)))
}

#simulations function: actually run the simulations
#Arguments: Args
#Returns: List with $fixed_results, $random_results
sims = function(args){
  checkArgs(args)
  random_results = NULL
  fixed_results = NULL
  for(s in 1:args$sims){
    Zomega = random(args$N, args$M, args$effect, args$h2)
    Xbeta = fixed(args$N, args$catCov, args$catBeta, args$binCov, args$binBeta, args$tdCov)
    if(args$dist[1]=='weibull'){
      data = weibull(args$N, Xbeta, Zomega, args$censor, args$dist[2], args$dist[3])
    } else if(args$dist[1]=='lognorm'){
      data = lognormal(args$N, Xbeta, Zomega, args$censor, args$dist[2], args$dist[3])
    }
    else{
      data = exponential(args$N, Xbeta, Zomega, args$censor, args$dist[2])
    }
    if(is.null(Zomega$GRM) & is.null(Xbeta$covars)) stop('There is no covariate or genetic effect on survival')
    if(is.null(Zomega$GRM) & !is.null(Xbeta$covars)){
      z = summary(coxph(Surv(data$y,data$death) ~ Xbeta$covars))
      fixed_results = rbind(fixed_results, z$coef[,1])
    } else if(!is.null(Zomega$GRM) & is.null(Xbeta$covars)) {
      z = coxmeg(outcome = cbind(data$y, data$death), corr = Zomega$GRM, type = 'dense', spd=F, solver = 'NM')
      random_results = c(random_results, z$tau)
    } else{
      z = coxmeg(outcome = cbind(data$y, data$death), corr = Zomega$GRM, type = 'dense', X = Xbeta$covars, spd=F, solver = 'NM')
      fixed_results = rbind(fixed_results, z$beta)
      random_results = c(random_results, z$tau)
    }
  } 
  if(!is.null(Xbeta$betas)){
    fixed_results = data.frame(fixed_results)
  }
  return(list("fixed_results"=fixed_results,"random_results"=random_results))
}

parser = parsing()
#All simulations:
# N=500, M=100, 100 sims, exponential distribution

#Sanity (NO Censoring)
#only genetic effect, random drawn multivariate normal h2 = 0.2, 0.8
onlyMVN20 = sims(parser$parse_args(c('-c', 0, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
onlyMVN80 = sims(parser$parse_args(c('-c', 0, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))

#only genetic effect, sum of fixed h2 = 0.2, 0.8
onlySUM20 = sims(parser$parse_args(c('-c', 0, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
onlySUM80 = sims(parser$parse_args(c('-c', 0, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))

#NO Censoring
#only binary fixed effect
#onlyBin1 = sims(parser$parse_args(c('-c', 0, '--binCov', 1, 0.5, '--binBeta', 0.7)))
#binary fixed effect + random drawn multivariate normal h2 = 0.2, 0.8
bin1MVN20 = sims(parser$parse_args(c('-c', 0, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
bin1MVN80 = sims(parser$parse_args(c('-c', 0, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))
#binary fixed effect + sum of fixed effect h2 = 0.2, 0.8
bin1SUM20 = sims(parser$parse_args(c('-c', 0, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
bin1SUM80 = sims(parser$parse_args(c('-c', 0, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))

#3 binary fixed effects
#onlyBin3 = sims(parser$parse_args(c('-c', 0, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4)))
#3 binary fixed effects + random drawn multivariate normal h2 = 0.2, 0.8
bin3MVN20 = sims(parser$parse_args(c('-c', 0, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
bin3MVN80 = sims(parser$parse_args(c('-c', 0, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))
#3 binary fixed effects + sum of fixed effect h2 = 0.2, 0.8
bin3SUM20 = sims(parser$parse_args(c('-c', 0, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
bin3SUM80 = sims(parser$parse_args(c('-c', 0, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))

#only categorical fixed effect
#onlyCat10 = sims(parser$parse_args(c('-c', 0, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77)))
#categorical fixed effect + random drawn multivariate normal h2 = 0.2, 0.8
cat10MVN20 = sims(parser$parse_args(c('-c', 0, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
cat10MVN80 = sims(parser$parse_args(c('-c', 0, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))
#categorical fixed effect + sum of fixed effect h2 = 0.2, 0.8
cat10SUM20 = sims(parser$parse_args(c('-c', 0, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
cat10SUM80 = sims(parser$parse_args(c('-c', 0, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))

#50% Censoring
#only binary fixed effect
#onlyBin1censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 1, 0.5, '--binBeta', 0.7)))
#binary fixed effect + random drawn multivariate normal h2 = 0.2, 0.8
bin1MVN20censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
bin1MVN80censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))
#binary fixed effect + SUM of fixed effect h2 = 0.2, 0.8
bin1SUM20censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
bin1SUM80censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 1, 0.5, '--binBeta', 0.7, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))

#3 binary fixed effects
#onlyBin3censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4)))
#3 binary fixed effects + random drawn multivariate normal h2 = 0.2, 0.8
bin3MVN20censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
bin3MVN80censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))
#3 binary fixed effects + sum of fixed effect h2 = 0.2, 0.8
bin3SUM20censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
bin3SUM80censor = sims(parser$parse_args(c('-c', 0.5, '--binCov', 3, 0.5, '--binBeta', 1.6, -0.7, 0.4, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))

#only categorical fixed effect
#onlyCat10censor = sims(parser$parse_args(c('-c', 0.5, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77)))
#categorical fixed effect + random drawn multivariate normal h2 = 0.2, 0.8
cat10MVN20censor = sims(parser$parse_args(c('-c', 0.5, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.2)))
cat10MVN80censor = sims(parser$parse_args(c('-c', 0.5, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'MVN', '-h2', 0.8)))
#categorical fixed effect + sum of fixed effect h2 = 0.2, 0.8
cat10SUM20censor = sims(parser$parse_args(c('-c', 0.5, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.2)))
cat10SUM80censor = sims(parser$parse_args(c('-c', 0.5, '--catCov', 10, '--catBeta', 0.39, 1.03, 0.69, -1.13, 0.35, 0.53, 0.61, -0.59, -0.77, '-d', 'weibull', 1.5, 0.5, '-g', 'SUM', '-h2', 0.8)))


final = data.frame(rep(c(0.2, 0.8), 12))
colnames(final) = "H2"
final["censoring"] = c(rep(0, 12), rep(0.5, 12))
final["geneticEffect"] = rep(c("MVN","MVN","SUM","SUM"),6)
final["fixedEffect"] = c(rep("Binary1",4),
                       rep("Binary3",4),
                       rep("Category10",4),
                       rep("Binary1",4),
                       rep("Binary3",4),
                       rep("Category10",4))
final["Mean"] = c(mean(bin1MVN20$random_results), mean(bin1MVN80$random_results), mean(bin1SUM20$random_results), mean(bin1SUM80$random_results),
                    mean(bin3MVN20$random_results), mean(bin3MVN80$random_results), mean(bin3SUM20$random_results), mean(bin3SUM80$random_results),
                    mean(cat10MVN20$random_results), mean(cat10MVN80$random_results), mean(cat10SUM20$random_results), mean(cat10SUM80$random_results),
                    mean(bin1MVN20censor$random_results), mean(bin1MVN80censor$random_results), mean(bin1SUM20censor$random_results), mean(bin1SUM80censor$random_results),
                    mean(bin3MVN20censor$random_results), mean(bin3MVN80censor$random_results), mean(bin3SUM20censor$random_results), mean(bin3SUM80censor$random_results),
                    mean(cat10MVN20censor$random_results), mean(cat10MVN80censor$random_results), mean(cat10SUM20censor$random_results), mean(cat10SUM80censor$random_results))
final["SE"] = c(sd(bin1MVN20$random_results), sd(bin1MVN80$random_results), sd(bin1SUM20$random_results), sd(bin1SUM80$random_results),
                    sd(bin3MVN20$random_results), sd(bin3MVN80$random_results), sd(bin3SUM20$random_results), sd(bin3SUM80$random_results),
                    sd(cat10MVN20$random_results), sd(cat10MVN80$random_results), sd(cat10SUM20$random_results), sd(cat10SUM80$random_results),
                    sd(bin1MVN20censor$random_results), sd(bin1MVN80censor$random_results), sd(bin1SUM20censor$random_results), sd(bin1SUM80censor$random_results),
                    sd(bin3MVN20censor$random_results), sd(bin3MVN80censor$random_results), sd(bin3SUM20censor$random_results), sd(bin3SUM80censor$random_results),
                    sd(cat10MVN20censor$random_results), sd(cat10MVN80censor$random_results), sd(cat10SUM20censor$random_results), sd(cat10SUM80censor$random_results))
final$SE = final$SE/sqrt(100)

ggplot(final[which(final$H2==0.2),],aes(x=censoring, y=Mean, fill = fixedEffect, color = geneticEffect, ymin = Mean-SE, ymax= Mean+SE)) +
  geom_bar(stat='identity', position = position_dodge(0.9)) + geom_errorbar(position = position_dodge(0.9)) + 

