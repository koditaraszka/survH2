
library(MASS)
library(coxmeg)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

source("~/Desktop/func_reml.R")
tol=1e-8

N = 500 # number of individuals
M = 100 # number of (causal) SNPs 
sims = 25 #as.numeric(args[3]) #number of simulations
herit = 50 #as.numeric(args[1])
h2 = herit/100
censoring = 0 #as.numeric(args[2])
cens = censoring/100
results = data.frame(rep(1:sims))
colnames(results) = "Sims"
results["H2"] = h2
results["censor"] = rep(cens, sims)
results["obsCensor"] = rep(NA, dim(results)[1])
results["QNORM"] = rep(NA, dim(results)[1])
#results["tauBOBYQA"] = rep(NA, dim(results)[1])
#results["tauBRENT"] = rep(NA, dim(results)[1])
results["tauNM"] = rep(NA, dim(results)[1])
results["varY"] = rep(NA, dim(results)[1])
results["varLogY"] = rep(NA, dim(results)[1])

for(s in 1:sims){
  MAF = runif(M, min=0.01, max=0.99)
  X = NULL
  for(m in 1:M){
    g = rbinom(N, 2, prob = MAF[m])
    X = cbind(X, g)
  }
  X = scale(X)
  GRM = X %*% t(X) / M
  beta = rnorm(M,0,1)
  xb = sqrt(h2)*scale(X %*% beta)
  #error = sqrt(1-h2)*scale(rnorm(N, 0, 1))
  #xb = mvrnorm(1, mu = rep(0,N), Sigma = h2*GRM)
  
  myrates = exp(xb)
  y_surv = rexp(N, rate = myrates)

  myrates2 = exp(xb-1)
  y_surv2 = rexp(N, rate = myrates2) #time to event
  if(cens==0){
    death = rep(1, N) # everyone has event
    log_y = log(y_surv)
    qnorm_y = qnorm(rank(y_surv, ties.method = "random")/(length(y_surv)+1))
    p_censor = length(which(death==0))/N
    varY = var(y_surv)
    varLogY = var(log_y)
    cheat = 
    #bobyqa_tau = coxmeg(cbind(y_surv, death), GRM, type = 'dense', spd=F)$tau
    nm_tau = coxmeg(cbind(y_surv, death), GRM, type = 'dense', spd=F, solver = 'NM')$tau
    #brent_tau = coxmeg(cbind(y_surv, death), GRM, type = 'dense', spd=F, solver = 'Brent')$tau
    qnorm = aiML(list(GRM), qnorm_y, c(0.5,0.5))$h2
  } else {
    cen <- rexp(N, rate = cens )
    y_cen <- pmin(y_surv, cen)
    death = as.numeric(y_surv <= cen)
    log_y = log(y_cen)
    qnorm_y = qnorm(rank(y_cen, ties.method = "random")/(length(y_cen)+1))
    p_censor = length(which(death==0))/N
    varY = var(y_cen)
    varLogY = var(log_y)
    #bobyqa_tau = coxmeg(cbind(y_surv, death), GRM, type = 'dense', spd=F)$tau
    nm_tau = coxmeg(cbind(y_surv, death), GRM, type = 'dense', spd=F, solver = 'NM')$tau
    #brent_tau = coxmeg(cbind(y_surv, death), GRM, type = 'dense', spd=F, solver = 'Brent')$tau
    qnorm = aiML(list(GRM), qnorm_y, c(0.5,0.5))$h2
  }
  results[which(results$Sims==s), "obsCensor"] = p_censor
  results[which(results$Sims==s), "QNORM"] = qnorm
  #results[which(results$Sims==s), "tauBOBYQA"] = bobyqa_tau
  #results[which(results$Sims==s), "tauBRENT"] = brent_tau
  results[which(results$Sims==s), "tauNM"] = nm_tau
  results[which(results$Sims==s), "varY"] = varY
  results[which(results$Sims==s), "varLogY"] = varLogY
  
  #write.table(results[which(results$Sims==s),], paste0("temp/temp_h2",herit,"_cens",censoring,".txt"), quote =F, col.names = F, row.names = F, sep = '\t', append=T)
}


final_coxmeg <- read.delim("~/Desktop/final_h2_coxmeg.txt", header=FALSE, stringsAsFactors = F)
colnames(final_coxmeg) = c("censor", paste0("obsCensor",c("Mean","SE")),
                           paste0("coxmeg",c("Mean","SE")),
                           paste0("coxmegPaper",c("Mean","SE")),
                           paste0("coxmegLog",c("Mean","SE")),
                           paste0("raw",c("Mean","SE")), 
                           paste0("log",c("Mean","SE")), 
                           paste0("qnorm",c("Mean","SE")), 
                           paste0("boxcox",c("Mean","SE")), 
                           paste0("probIntegral",c("Mean","SE")),
                           paste0("sanity",c("Mean","SE")), "H2")

final_coxmeg$obsCensorMean = round(final_coxmeg$obsCensorMean, 3)
final_coxmeg$obsCensorSE = round(final_coxmeg$obsCensorSE, 4)
paste("Coxmeg:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$coxmegMean + final_coxmeg$obsCensorMean + final_coxmeg$coxmegMean*final_coxmeg$obsCensorMean))$r.squared,3))
paste("Coxmeg Paper:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$coxmegPaperMean + final_coxmeg$obsCensorMean + final_coxmeg$coxmegPaperMean*final_coxmeg$obsCensorMean))$r.squared,3))
paste("Qnorm:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$qnormMean + final_coxmeg$obsCensorMean + final_coxmeg$qnormMean*final_coxmeg$obsCensorMean))$r.squared,3))
paste("Raw:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$rawMean + final_coxmeg$obsCensorMean + final_coxmeg$rawMean*final_coxmeg$obsCensorMean))$r.squared,3))

paste("Log:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$logMean + final_coxmeg$obsCensorMean + final_coxmeg$logMean*final_coxmeg$obsCensorMean))$r.squared,3))
paste("BoxCox:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$boxcoxMean + final_coxmeg$obsCensorMean + final_coxmeg$boxcoxMean*final_coxmeg$obsCensorMean))$r.squared,3))
paste("Prob Integral:", round(summary(lm(final_coxmeg$H2 ~ final_coxmeg$probIntegralMean + final_coxmeg$obsCensorMean + final_coxmeg$probIntegralMean*final_coxmeg$obsCensorMean))$r.squared,3))

final_figdata = data.frame(rep(final_coxmeg$H2, 5))
colnames(final_figdata) = "True Heritability"
final_figdata$`True Heritability` = factor(final_figdata$`True Heritability`, levels = unique(final_coxmeg$H2))
final_figdata["Censoring"] = rep(final_coxmeg$censor, 5)
final_figdata["ObsCensoring"] = rep(final_coxmeg$obsCensorMean, 5)
final_figdata$Censoring = factor(final_figdata$Censoring, levels = unique(final_coxmeg$censor))
final_figdata["Transformation"] = c(rep("None",dim(final_coxmeg)[1]), rep("Sanity",dim(final_coxmeg)[1]), rep("Qnorm",dim(final_coxmeg)[1]), rep("COXMEG",dim(final_coxmeg)[1]), rep("COXMEG_Paper",dim(final_coxmeg)[1]))
final_figdata["Estimated Heritability"] = c(final_coxmeg$rawMean, final_coxmeg$sanityMean, final_coxmeg$qnormMean, final_coxmeg$coxmegMean, final_coxmeg$coxmegPaperMean)
final_figdata["SE"] = c(final_coxmeg$rawSE, final_coxmeg$sanitySE, final_coxmeg$qnormSE, final_coxmeg$coxmegSE, final_coxmeg$coxmegPaperSE)

ggplot(final_figdata, aes(x=`True Heritability`, y=`Estimated Heritability`, fill = Transformation, ymin = `Estimated Heritability`-SE, ymax = `Estimated Heritability`+SE)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + facet_wrap(~Censoring) + geom_errorbar(position = position_dodge(0.9)) + theme_minimal()

ggplot(final_figdata[which(final_figdata$Transformation=="COXMEG_Paper"),], aes(x=`True Heritability`, y=`Estimated Heritability`, fill = factor(Censoring), ymin = `Estimated Heritability`-SE, ymax = `Estimated Heritability`+SE)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + facet_wrap(~Transformation) + geom_errorbar(position = position_dodge(0.9)) + theme_minimal()


fig_data <- read.delim("~/Desktop/fig_data.txt", header=F, stringsAsFactors = F)
colnames(fig_data) = c("Sims", "H2", "censor", "obsCensor", "QNORM", "tauBOBYQA", 
                       "tauBRENT", "tauNM", "varY", "varLogY")

fig_data = fig_data %>% group_by(H2, censor) %>% summarise(
  meanObsCensor = mean(obsCensor),
  meanQNORM = mean(QNORM),
  meanTauNM = mean(tauNM),
  meanTauVarY = mean(tauNM/varY),
  meanTauVarLogY = mean(tauNM/varLogY),
  meanPaperTrueCensor = mean(tauNM/(tauNM + 1/(1 - censor))),
  meanPaperObsCensor = mean(tauNM/(tauNM + 1/(1 - obsCensor))),
  seObsCensor = sd(obsCensor)/sqrt(50),
  seQNORM = sd(QNORM)/sqrt(50),
  seTauNM = sd(tauNM)/sqrt(50),
  seTauVarY = sd(tauNM/varY)/sqrt(50),
  seTauVarLogY = sd(tauNM/varLogY)/sqrt(50),
  sePaperTrueCensor = sd(tauNM/(tauNM + 1/(1 - censor)))/sqrt(50),
  sePaperObsCensor = sd(tauNM/(tauNM + 1/(1 - obsCensor)))/sqrt(50))

final_figdata = data.frame(rep(fig_data$H2, 4))
colnames(final_figdata) = "True Heritability"
final_figdata$`True Heritability` = factor(final_figdata$`True Heritability`, levels = unique(fig_data$H2))
final_figdata["Censoring"] = rep(fig_data$censor, 4)
final_figdata["ObsCensoring"] = rep(fig_data$meanObsCensor, 4)
final_figdata$Censoring = factor(final_figdata$Censoring, levels = unique(fig_data$censor))
final_figdata["Transformation"] = c(rep("QNORM",dim(fig_data)[1]), rep("Tau",dim(fig_data)[1]), 
                                    rep("PaperTrueCensor",dim(fig_data)[1]), rep("PaperObsCensor",dim(fig_data)[1]))
final_figdata["Estimated Heritability"] = c(fig_data$meanQNORM, fig_data$meanTauNM, fig_data$meanPaperTrueCensor, fig_data$meanPaperObsCensor) 
final_figdata["SE"] = c(fig_data$seQNORM, fig_data$seTauNM, fig_data$sePaperTrueCensor, fig_data$sePaperObsCensor)

ggplot(final_figdata, aes(x=`True Heritability`, y=`Estimated Heritability`, fill = Transformation, ymin = `Estimated Heritability`-SE, ymax = `Estimated Heritability`+SE)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + facet_wrap(~Censoring) + geom_errorbar(position = position_dodge(0.9)) + theme_minimal()

ggplot(final_figdata[which(final_figdata$Transformation=="COXMEG_Paper"),], aes(x=`True Heritability`, y=`Estimated Heritability`, fill = factor(Censoring), ymin = `Estimated Heritability`-SE, ymax = `Estimated Heritability`+SE)) + 
  geom_bar(stat = 'identity', position = position_dodge(0.9)) + facet_wrap(~Transformation) + geom_errorbar(position = position_dodge(0.9)) + theme_minimal()


final_figdata = data.frame(rep(finalfixed_data$h2, 4))
colnames(final_figdata) = "True Heritability"
final_figdata$`True Heritability` = factor(final_figdata$`True Heritability`, levels = unique(finalfixed_data$h2))
final_figdata["Censoring"] = rep(finalfixed_data$censor, 4)
final_figdata["ObsCensoring"] = rep(finalfixed_data$meanObsCensor, 4)
final_figdata$Censoring = factor(final_figdata$Censoring, levels = unique(finalfixed_data$censor))
final_figdata["Transformation"] = c(rep("QNORM",dim(finalfixed_data)[1]), rep("Tau",dim(finalfixed_data)[1]), 
                                    rep("PaperTrueCensor",dim(finalfixed_data)[1]), rep("PaperObsCensor",dim(finalfixed_data)[1]))
final_figdata["Estimated Heritability"] = c(finalfixed_data$meanQNORM, finalfixed_data$meanTauNM, finalfixed_data$meanPaperTrueCensor, finalfixed_data$meanPaperObsCensor) 
final_figdata["SE"] = c(finalfixed_data$seQNORM, finalfixed_data$seTauNM, finalfixed_data$sePaperTrueCensor, finalfixed_data$sePaperObsCensor)



finalrandom_data = final_coxmeg_random %>% group_by(h2, censor) %>% summarise(
  meanObsCensor = mean(obsCensor),
  meanQNORM = mean(QNORM),
  meanTauNM = mean(tauNM),
  meanTauVarY = mean(tauNM/varY),
  meanTauVarLogY = mean(tauNM/varLogY),
  meanPaperTrueCensor = mean(tauNM/(tauNM + 1/(1 - censor))),
  meanPaperObsCensor = mean(tauNM/(tauNM + 1/(1 - obsCensor))),
  seObsCensor = sd(obsCensor)/sqrt(50),
  seQNORM = sd(QNORM)/sqrt(50),
  seTauNM = sd(tauNM)/sqrt(50),
  seTauVarY = sd(tauNM/varY)/sqrt(50),
  seTauVarLogY = sd(tauNM/varLogY)/sqrt(50),
  sePaperTrueCensor = sd(tauNM/(tauNM + 1/(1 - censor)))/sqrt(50),
  sePaperObsCensor = sd(tauNM/(tauNM + 1/(1 - obsCensor)))/sqrt(50))


finalfixed_data = final_coxmeg_binomial %>% group_by(h2, censor) %>% summarise(
  meanObsCensor = mean(obsCensor),
  meanQNORM = mean(QNORM),
  meanTauNM = mean(tauNM),
  meanTauVarY = mean(tauNM/varY),
  meanTauVarLogY = mean(tauNM/varLogY),
  meanPaperTrueCensor = mean(tauNM/(tauNM + 1/(1 - censor))),
  meanPaperObsCensor = mean(tauNM/(tauNM + 1/(1 - obsCensor))),
  seObsCensor = sd(obsCensor)/sqrt(50),
  seQNORM = sd(QNORM)/sqrt(50),
  seTauNM = sd(tauNM)/sqrt(50),
  seTauVarY = sd(tauNM/varY)/sqrt(50),
  seTauVarLogY = sd(tauNM/varLogY)/sqrt(50),
  sePaperTrueCensor = sd(tauNM/(tauNM + 1/(1 - censor)))/sqrt(50),
  sePaperObsCensor = sd(tauNM/(tauNM + 1/(1 - obsCensor)))/sqrt(50))
