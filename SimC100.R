library(pracma) 
dyn.load("src.dll")
source("KC_correction.R")
source("Etai2Deltai1.R")
source("GenBinaryYGamma.R")
source("get.GH1.R")

#install.packages('devtools')
#library(devtools)
#install_github("mercaldo/MMLB", force=TRUE)
library("geepack")
library("geesmv")
library(MMLB)

set.seed(1)
R = 1000 # replication times
p = 0.5 # proportion of treatment group
N = 100 # number of cluster
n = 100 # cluster size
beta=c(-0.8, 1.2, 1.5, 0.4)

GEE_exch_coef <- NULL # store coefficient results for GEE working eachangeable
MMMcoef <- NULL # store coefficient results for MMM
GEE_exch_sd <- NULL #store estimated coefficient standard deviation results for GEE working eachangeable
GEE_cor_sd <- NULL # store MD corrected coefficient standard deviation estimation for GEE working exchangeable
GEE_cor_kc <- NULL # store KC corrected coefficient standard deviation estimation for GEE working exchangeable
GEE_cor_fg <- NULL # store FG corrected coefficient standard deviation estimation for GEE working exchangeable
MMMsd <- NULL # store estimated coefficient standard deviation results for MMM

for(i in 1:R){
nclust = rep(n,N)
id     = rep(seq(N), nclust)
trt    =rep(0,N)
trt[1:(N*p)]=1
Xe     = rep(trt, nclust) # binary exposure
time   = as.numeric(sapply( nclust, function(ZZ) {rnorm(ZZ)} ))
data   = data.frame(id, time, Xe)
data   = data[order(data$id),]

newdata1 = GenBinaryYGamma(mean.formula=~time*Xe, lv.formula=~1, 
          beta=beta, sigma=0.5, lambda=1, id=id, data=data, q=120, 
          Yname = "binY")$data
wts <- rep(1, nrow(newdata1))
newdata1 <- cbind(newdata1, wts)
tryCatch({
gee.fit <-  geeglm(binY~time+Xe+time:Xe, id = id, data = newdata1, family = binomial, corstr = "exchangeable", scale.fix=T)
GEE_exch_coef <- rbind(GEE_exch_coef, c(i, coef(gee.fit)))
GEE_exch_sd <- rbind(GEE_exch_sd,c(i, summary(gee.fit)$coefficients[,2]))
}, warning=function(w){print(paste("GEE_ind",i,w))}, error=function(e){print(paste("GEE_ind",i,e))})

tryCatch({
md.exch <- GEE.var.md(binY~time+Xe+time:Xe,id=id,family="binomial", data = newdata1,corstr="exchangeable") 
GEE_cor_sd <- rbind(GEE_cor_sd,c(i, sqrt(md.exch$cov.beta)))
kc.exch <- GEE.var.kc3(data=newdata1, model=as.formula("binY~time+Xe+time:Xe"), betahat=coefficients(gee.fit), phi=1, 
                  alpha = as.numeric(summary(gee.fit)$corr[1]))
GEE_cor_kc <- rbind(GEE_cor_kc, c(i, kc.exch))
fg.exch <- GEE.var.fg(binY~time+Xe+time:Xe,id=id,family="binomial", data = newdata1,corstr="exchangeable") 
GEE_cor_fg <- rbind(GEE_cor_fg,c(i, sqrt(fg.exch$cov.beta)))
}, warning=function(w){print(paste("GEE_ex",i,w))}, error=function(e){print(paste("GEE_ex",i,e))})

tryCatch({
library("MMLB")
mod_mtlv1 = mm(binY~time*Xe,lv.formula=~1, data=newdata1,id=id, q=120)
MMMcoef <- rbind(MMMcoef, c(i, summary(mod_mtlv1)[[5]][,1]))
MMMsd <- rbind(MMMsd, c(i, summary(mod_mtlv1)[[5]][,2]))
}, warning=function(w){print(paste("MMM",i,w))}, error=function(e){print(paste("MMM",i,e))})
}

write.csv(GEE_exch_coef, paste0("GEE_exchangeable_coef_K", N, "n",n, ".csv"))
write.csv(GEE_exch_sd, paste0("GEE_exchangeable_sd_K", N, "n",n, ".csv"))
write.csv(GEE_cor_sd, paste0("GEE_independence_sd_K", N, "n",n, ".csv"))
write.csv(GEE_cor_kc, paste0("GEE_kc_sd_K", N, "n",n, "1.csv"))
write.csv(GEE_cor_fg, paste0("GEE_fg_sd_K", N, "n",n, "1.csv"))
write.csv(MMMcoef, paste0("MMMcoef_K", N, "n",n, ".csv"))
write.csv(MMMsd, paste0("MMMsd_K", N, "n",n, ".csv"))

