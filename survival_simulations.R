library(MASS)
library(coxed)
library(coxme)
library(coxmeg)
library(dplyr)
library(ggplot2)

#args = commandArgs(trailingOnly=TRUE)
rm(list=ls())
gc()

source("~/Desktop/func_reml.R")
tol=1e-8

N = 500 # number of individuals
M = 100 # number of (causal) SNPs 
sims = 50 #number of simulations
heritability = from 0.06 to 0.5, step of 0.02#seq(from=100, to=500, by=25)/1000
censoring = from 0.0 to 0.15, ste of 0.01#seq(from=2, to=15, by =1)/100

final = NULL
for(h2 in heritability){
  results = data.frame(rep(1:sims,(length(censoring)+1))) #data.frame(rep(1:sims,2*(length(censoring)+1)))
  colnames(results) = "Sims"
  results["Method"] = c(rep("REML",sims*(length(censoring)+1))) #rep("COXME",sims*(length(censoring)+1)))
  censor = rep(0, sims)
  for(cens in censoring){
    censor = c(censor, rep(cens, sims))
  }
  results["censor"] = censor #c(censor, rep(censor,2)
  
  results["raw"] = rep(NA, dim(results)[1])
  results["log"] = rep(NA, dim(results)[1])
  results["qnorm"] = rep(NA, dim(results)[1])
  results["sanity"] = rep(NA, dim(results)[1])
  results["boxcox"] = rep(NA, dim(results)[1])
  results["probIntegral"] = rep(NA, dim(results)[1])
  for(s in 1:sims){
    X = mvrnorm(N, mu = rep(0,M), Sigma = diag(M)) # change the binomial MAF from 0.01, 0.5 
    GRM = X %*% t(X) / M
    beta = rnorm(M,0,1)
    xb = sqrt(h2)*scale(X %*% beta)
  
    myrates = exp(xb)
    y_surv = rexp(N, rate = myrates) #time to event
    log_y = log(y_surv)
    qnorm_y = qnorm(rank(y_surv, ties.method = "random")/(length(y_surv)+1))
    death = rep(1, N) # everyone has event
    
    error = sqrt(1-h2)*scale(rnorm(N, 0, 1))
    y_sanity = xb + error
    
    bc = boxcox(lm(y_surv ~ 1), plotit =F)
    lambda <- bc$x[which.max(bc$y)]
    boxcox_y <- (y_surv ^ lambda - 1) / lambda
    u = 1 - exp(-(1/mean(y_surv))*y_surv)
    u[which(u==1)] = 1-tol
    u[which(u==0)] = tol
    z = qnorm(u)
    results[which(results$Sims==s & results$censor==0), "sanity"] = aiML(list(GRM), y_sanity, c(0.5,0.5))$h2
    results[which(results$Sims==s & results$censor==0), "raw"] = aiML(list(GRM), y_surv, c(0.5,0.5))$h2
    results[which(results$Sims==s &results$censor==0), "log"] = aiML(list(GRM), log_y, c(0.5,0.5))$h2
    results[which(results$Sims==s &results$censor==0), "qnorm"] = aiML(list(GRM), qnorm_y, c(0.5,0.5))$h2
    results[which(results$Sims==s & results$censor==0), "boxcox"] = aiML(list(GRM), boxcox_y, c(0.5,0.5))$h2
    results[which(results$Sims==s & results$censor==0), "probIntegral"] = aiML(list(GRM), z, c(0.5,0.5))$h2
    
    for(cens in censoring){
      cen <- rexp(N, rate = cens )
      y_cen <- pmin(y_surv, cen)
      death = as.numeric(y_surv <= cen)
      
      
      y_cens = y_surv[which(y_surv <= cen)]
      #can you check the fraction of that are censored are actually being censored
      #15% set to be censored but only 1% censored?
      #can we look at the simulations based on the observed censoring rate not the true censoring
      log_y = log(y_cen)
      qnorm_y = qnorm(rank(y_cen, ties.method = "random")/(length(y_cen)+1))
      bc = boxcox(lm(y_cen ~ 1), plotit=F)
      lambda <- bc$x[which.max(bc$y)]
      boxcox_y <- (y_cen ^ lambda - 1) / lambda
      u = 1 - exp(-(1/mean(y_cen))*y_cen)
      u[which(u>(1-tol))] = 1-tol
      u[which(u<tol)] = tol
      z = qnorm(u)
      
      results[which(results$Sims==s & results$censor==cens), "raw"] = aiML(list(GRM), y_cen, c(0.5,0.5))$h2
      results[which(results$Sims==s &results$censor==cens), "log"] = aiML(list(GRM), log_y, c(0.5,0.5))$h2
      results[which(results$Sims==s &results$censor==cens), "qnorm"] = aiML(list(GRM), qnorm_y, c(0.5,0.5))$h2
      results[which(results$Sims==s &results$censor==cens), "boxcox"] = aiML(list(GRM), boxcox_y, c(0.5,0.5))$h2
      results[which(results$Sims==s & results$censor==cens), "probIntegral"] = aiML(list(GRM), z, c(0.5,0.5))$h2
      
    }
  }
  temp = data.frame(results %>% group_by(censor) %>% summarize(rawMean = mean(raw), rawSE = sd(raw)/sqrt(N), 
                                                               logMean = mean(log), logSE = sd(log)/sqrt(N), 
                                                               qnormMean = mean(qnorm), qnormSE = sd(qnorm)/sqrt(N),
                                                               boxcoxMean = mean(boxcox, na.rm=T), boxcoxSE = sd(boxcox, na.rm=T)/sqrt(N),
                                                               probIntegralMean = mean(probIntegral, na.rm=T), probIntegralSE = sd(probIntegral, na.rm=T)/sqrt(N),
                                                               sanityMean = mean(sanity, na.rm=T), sanitySE = sd(sanity, na.rm=T)/sqrt(N)))
  temp["H2"] = h2
  final = rbind(final, temp)
}

fig_data = data.frame(rep(final_downsample$H2, 6))
colnames(fig_data) = "True Heritability"
fig_data$`True Heritability` = factor(fig_data$`True Heritability`, levels = heritability)
fig_data["Censoring"] = rep(final_downsample$censor, 6)
fig_data$Censoring = factor(fig_data$Censoring, levels = c(0, censoring))
fig_data["Transformation"] = c(rep("None",dim(final_downsample)[1]),rep("Log",dim(final_downsample)[1]),rep("Qnorm",dim(final_downsample)[1]),
                               rep("BoxCox",dim(final_downsample)[1]),rep("ProbIntegral",dim(final_downsample)[1]),rep("Sanity",dim(final_downsample)[1]))
fig_data["Estimated Heritability"] = c(final_downsample$rawMean, final_downsample$logMean, final_downsample$qnormMean,
                                       final_downsample$boxcoxMean, final_downsample$probIntegralMean, final_downsample$sanityMean)
fig_data["SE"] = c(final_downsample$rawSE, final_downsample$logSE, final_downsample$qnormSE,
                                       final_downsample$boxcoxSE, final_downsample$probIntegralSE, final_downsample$sanitySE)

ggplot(fig_data, aes(x=`True Heritability`, y=`Estimated Heritability`, fill = Transformation, ymin = `Estimated Heritability`-SE, ymax = `Estimated Heritability`+SE)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + facet_wrap(~Censoring) + geom_errorbar(position = position_dodge(0.9)) + theme_minimal()

ggplot(fig_data[-which(fig_data$Transformation=="Sanity"),], aes(x=`True Heritability`, y=`Estimated Heritability`, fill = factor(Censoring), ymin = `Estimated Heritability`-SE, ymax = `Estimated Heritability`+SE)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + facet_wrap(~Transformation) + geom_errorbar(position = position_dodge(0.9)) + theme_minimal()


final_downsample <- read.delim("~/Desktop/final_h2_downsample.txt", header=FALSE, stringsAsFactors = F)
colnames(final_downsample) = c("censor",paste0("raw",c("Mean","SE")), 
                                       paste0("log",c("Mean","SE")), 
                                       paste0("qnorm",c("Mean","SE")), 
                                       paste0("boxcox",c("Mean","SE")), 
                                       paste0("probIntegral",c("Mean","SE")),
                                       paste0("sanity",c("Mean","SE")), "H2")

paste("Raw:", round(summary(lm(final_downsample$H2 ~ final_downsample$rawMean + final_downsample$censor + final_downsample$rawMean*final_downsample$censor))$r.squared,3))
paste("Log:", round(summary(lm(final_downsample$H2 ~ final_downsample$logMean + final_downsample$censor + final_downsample$logMean*final_downsample$censor))$r.squared,3))
paste("Qnorm:", round(summary(lm(final_downsample$H2 ~ final_downsample$qnormMean + final_downsample$censor + final_downsample$qnormMean*final_downsample$censor))$r.squared,3))
paste("BoxCox:", round(summary(lm(final_downsample$H2 ~ final_downsample$boxcoxMean + final_downsample$censor + final_downsample$boxcoxMean*final_downsample$censor))$r.squared,3))
paste("Prob Integral:", round(summary(lm(final_downsample$H2 ~ final_downsample$probIntegralMean + final_downsample$censor + final_downsample$probIntegralMean*final_downsample$censor))$r.squared,3))


final_censor <- read.delim("~/Desktop/final_h2_censor.txt", header=FALSE, stringsAsFactors = F)
colnames(final_censor) = c("censor",paste0("raw",c("Mean","SE")), 
                               paste0("log",c("Mean","SE")), 
                               paste0("qnorm",c("Mean","SE")), 
                               paste0("boxcox",c("Mean","SE")), 
                               paste0("probIntegral",c("Mean","SE")),
                               paste0("sanity",c("Mean","SE")), "H2")

paste("Raw:", round(summary(lm(final_censor$H2 ~ final_censor$rawMean + final_censor$censor + final_censor$rawMean*final_censor$censor))$r.squared,3))
paste("Log:", round(summary(lm(final_censor$H2 ~ final_censor$logMean + final_censor$censor + final_censor$logMean*final_censor$censor))$r.squared,3))
paste("Qnorm:", round(summary(lm(final_censor$H2 ~ final_censor$qnormMean + final_censor$censor + final_censor$qnormMean*final_censor$censor))$r.squared,3))
paste("BoxCox:", round(summary(lm(final_censor$H2 ~ final_censor$boxcoxMean + final_censor$censor + final_censor$boxcoxMean*final_censor$censor))$r.squared,3))
paste("Prob Integral:", round(summary(lm(final_censor$H2 ~ final_censor$probIntegralMean + final_censor$censor + final_censor$probIntegralMean*final_censor$censor))$r.squared,3))
