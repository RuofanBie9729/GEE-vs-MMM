#install.packages('devtools')
#library(devtools)
#install_github("mercaldo/MMLB", force=TRUE)
library("MMLB")
library("geepack")
library("geesmv")
library(aod)
source("KC_correction.R")

set.seed(1)
R = 1000 # replication times
p = 0.5 # proportion of treatment group
N = 100 # number of cluster
n = 100 # cluster size
beta=c(-0.8, 1.2, 0)

GEE_exch <- NULL # store wald test statistics and p values for GEE working exchangeable 
MMM <- NULL # store wald test statistics and p values for MMM
GEE_cor <- NULL # store MD corrected wald test statistics and p values for GEE working exchangeable
MMM_cor <- NULL # store DF corrected wald test statistics and p values for MMM
MMM_lrt <- NULL # store likelihood ratio test results for MMM
GEE_cor_kc <- NULL # store KC corrected wald test statistics and p values for GEE working exchangeable
GEE_cor_fg <- NULL # store FG corrected wald test statistics and p values for GEE working exchangeable

for(i in 1:R){
nclust = rep(n,N)
id     = rep(seq(N), nclust)
trt    =rep(0,N)
trt[1:(N*p)]=1
Xe     = rep(trt, nclust) # binary exposure
time   = as.numeric(sapply( nclust, function(ZZ) {rnorm(ZZ)} ))
data   = data.frame(id, time, Xe)
data   = data[order(data$id, data$time),]

newdata1 = GenBinaryY(mean.formula=~time+Xe, lv.formula=~1, 
          beta=beta, sigma=0.5, id=id, data=data, q=120, 
          Yname = "binY")
wts <- rep(1, nrow(newdata1))
newdata1 <- cbind(newdata1, wts)
tryCatch({ ###GEE wald
gee.fit <-  geeglm(binY~time+Xe, id = id, data = newdata1, family = binomial, corstr = "exchangeable")
GEE_exch <- rbind(GEE_exch, c(i, summary(gee.fit)[[6]][,3] , summary(gee.fit)[[6]][,4]))
}, warning=function(w){print(paste("GEE_ind",i,w))}, error=function(e){print(paste("GEE_ind",i,e))})

tryCatch({ ###GEE correction wald
md.exch <- GEE.var.md(binY~time+Xe,id=id,family=binomial, data = newdata1,corstr="exchangeable") 
wald_stat <- (summary(gee.fit)[[6]][,1]/sqrt(md.exch$cov.beta))^2
p_value <- 1-pchisq(wald_stat, 1)
GEE_cor <- rbind(GEE_cor,c(i, wald_stat, p_value))
kc.exch <- GEE.var.kc3(data=newdata1, model=as.formula("binY~time+Xe"), betahat=coefficients(gee.fit), phi=1, 
                  alpha = as.numeric(summary(gee.fit)$corr[1]))
wald_stat_kc <- (summary(gee.fit)[[6]][,1]/kc.exch)^2
p_value_kc <- 1-pchisq(wald_stat_kc, 1)
GEE_cor_kc <- rbind(GEE_cor_kc, c(i, wald_stat_kc, p_value_kc))
fg.exch <- GEE.var.fg(binY~time+Xe,id=id,family="binomial", data = newdata1,corstr="exchangeable") 
wald_stat_fg <- (summary(gee.fit)[[6]][,1]/sqrt(fg.exch$cov.beta))^2
p_value_fg <- 1-pchisq(wald_stat_fg, 1)
GEE_cor_fg <- rbind(GEE_cor_fg,c(i, wald_stat_fg, p_value_fg))
}, warning=function(w){print(paste("GEE_ex",i,w))}, error=function(e){print(paste("GEE_ex",i,e))})

tryCatch({ ###MMM wald
mod_mtlv1 = mm(binY~time+Xe,lv.formula=~1, data=newdata1,id=id, q=120)
MMM <- rbind(MMM, c(i, summary(mod_mtlv1)[[5]][,3], summary(mod_mtlv1)[[5]][,4]))
}, warning=function(w){print(paste("MMM",i,w))}, error=function(e){print(paste("MMM",i,e))})

tryCatch({ ###MMM correction wald
wald_stat <- summary(mod_mtlv1)[[5]][,3]
p_value <- 1-pf(wald_stat, 1, N-2)
MMM_cor <- rbind(MMM_cor, c(i, wald_stat, p_value))
}, warning=function(w){print(paste("MMM",i,w))}, error=function(e){print(paste("MMM",i,e))})

tryCatch({ ###MMM LRT
mod_mtlv2 = mm(binY~time,lv.formula=~1, data=newdata1,id=id, q=120)
LRT_stat  <- -2*(mod_mtlv2$logLik-mod_mtlv1$logLik)
p_value <- 1-pchisq(LRT_stat, 1)
MMM_lrt <- rbind(MMM_lrt, c(i, LRT_stat, p_value))
}, warning=function(w){print(paste("MMM",i,w))}, error=function(e){print(paste("MMM",i,e))})

}

write.csv(GEE_exch, paste0("GEE_exchangeable_K", N, "n",n, "H0", ".csv"))
write.csv(GEE_cor, paste0("GEE_correction_K", N, "n",n, "H0", ".csv"))
write.csv(GEE_cor_kc, paste0("GEE_kc_sd_K", N, "n",n, "H0", ".csv"))
write.csv(GEE_cor_fg, paste0("GEE_fg_sd_K", N, "n",n, "H0", ".csv"))
write.csv(MMM, paste0("MMM_K", N, "n",n, "H0", ".csv"))
write.csv(MMM_cor, paste0("MMM_correction_K", N, "n",n, "H0", ".csv"))
write.csv(MMM_lrt, paste0("MMM_lrt_K", N, "n",n, "H0", ".csv"))

