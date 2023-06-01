# This is sample code for the approach suggested in Lo, S.M.S., Wilke, R.A. and Emura, T. (2023): A semiparametric model for the cause-specific hazard under risk proportionality, Working Paper.
# In particular, it produces the results of the application in Section 5 of the paper.

# It requires LFS.csv, the data file.
# The path for reading the data and for saving the figures needs to be set by the user.


rm(list = ls())
library(copula)
library(survival)
library(pastecs)
library(cmprsk)
library(timereg)
library(MALDIquant)
library(plyr)
library(dplyr)
library(labelled)
library(rms)
library(survival)
library("survminer")
library(stats4)
library(bbmle)
library(latex2exp)

############ Define several functions ##############

invpsigumbel <- function(s, tau){  # inverse copula generator of gumbel
  if (tau==0){
    u <- exp(-s)
  }  
  if (tau!=0){
    theta <- iTau(archmCopula("gumbel"), tau)
    u <- exp(-s^(1/theta))
  }
  return(u)
}

psigumbel<-function(s,tau){ # copula generator of gumbel
  if (tau==0){
    
    u<- -log(s)
  }
  if (tau!=0){
    theta <- iTau(archmCopula("gumbel"), tau)
    
    u<- (-log(s))^(theta)
  }
  return(u)
}

################# 

data0 <- read.csv("PATH/LFS.csv") ## input data file

x<-data0$wdur # input the name of duration variable 
z<-cbind(rep(1,length(x)),data0$jobno2, data0$twdur) # input the name of covariate
del = 2-data0$single  # input the name of risk indicator risk 1 = 1, risk 2 = 0  

#################
y<-data.frame(x,del,z)
summary(y)
nrow(y)

#Choose parametric model to be fitted
model = "weib" ## choose the model "expo" "weib" 


##################################################
# Data reorganisation for partial likelihood (28)
# replicate the data: need a y data frame y<-data.frame(x,del,z),
# x is the duration
# z includes covaraite without constant
# risk indicator del = 1 (risk 1) or 2 (risk 2) , no censoring

y1<-y[which(y[,2]==1),] 
y2<-y[which(y[,2]==2),] 
summary(y1)
nrow(y1)
summary(y2)
nrow(y2)

z11<-cbind(y1[,3:ncol(y1)],matrix(0,nrow=nrow(y1), ncol =ncol(y1)-2))
z12<-cbind(matrix(0,nrow=nrow(y1), ncol =ncol(y1)-2),y1[,3:ncol(y1)])
z21<-cbind(y2[,3:ncol(y2)],matrix(0,nrow=nrow(y2), ncol =ncol(y2)-2))
z22<-cbind(matrix(0,nrow=nrow(y2), ncol =ncol(y2)-2),y2[,3:ncol(y2)])

y11<-cbind(y1[,1],rep(1,nrow(y1)),z11)
y12<-cbind(y1[,1],rep(0,nrow(y1)),z12)
y21<-cbind(y2[,1],rep(0,nrow(y2)),z21)
y22<-cbind(y2[,1],rep(1,nrow(y2)),z22)

colnames(y11)<-c('X','del','z0','z1')
colnames(y12)<-c('X','del','z0','z1')
colnames(y21)<-c('X','del','z0','z1')
colnames(y22)<-c('X','del','z0','z1')

Y<- rbind(y11,y12,y21,y22)
colnames(Y) <- c('x','del','a1','b11','b12', 'a2','b21','b22')

#### this uses the coxph to estimate the partial likelihood estimates

res<-coxph(Surv(Y[,1],Y[,2]) ~  Y$b11 + Y$b12 +Y$a2 +Y$b21+Y$b22) # add covariates if necessary
res
#c(b1,a2,b2) # it is to double check with the true data


Y2<- cbind(Y[,1],Y[,2],Y[,3],Y[,4],Y[5]/52,Y[,6],Y[7],Y[8]/52)

colnames(Y2) <- c('x','del','a1','b11','b12', 'a2','b21','b22')

res2<-coxph(Surv(Y2[,1],Y2[,2]) ~  Y2$b11 + Y2$b12 +Y2$a2 +Y2$b21+Y2$b22) # add covariates if necessary
res2

test<-cox.zph(res) # it is the schoenfeld test
print(test) # Ho: linar, p>0.05 means not rejected.

#### estimate alpha 1 - the constant

bcumh <- basehaz(res,centered=FALSE) # it can be use to check the baseline cum. hazard 
fit<-lm(formula = bcumh[which(bcumh[,2]<600),1] ~ -1 + bcumh[which(bcumh[,2]<600),2] ) # regress integrated varsigma(t) on t, nn constant
summary(fit)

bh<-fit$coefficients[[1]] # estimate the baseline cumulative hazard

plot(bcumh[which(bcumh[,2]<600),2],bcumh[which(bcumh[,2]<600),1], type="l", lty=1, ylab="cumulative hazard", xlab="Job duration") 
abline(a=0, b=bh, lty=2)
legend("bottom",legend=c("estimated baseline cumulative hazard", "fit for exponential hazard"), lty=1:2)

ea1<-fit$coefficients[[1]]  #### estimated alpha 1
ea1

eb11<-res2$coefficients[[1]] # estiamted beta11
eb12<-res2$coefficients[[2]]  # estiamted beta12
ea2<-exp(res2$coefficients[[3]])   # estiamted alpha2
eb21<-res2$coefficients[[4]]   # estiamted beta21
eb22<-res2$coefficients[[5]]   # estiamted beta22


#############################################

#data preparation of original sample
t<-data0$wdur # input the name of duration variable 
z1<-data0$jobno2 #covariates
z2<- data0$twdur/52 
delw = 2-data0$single  # input the name of risk indicator risk 1 = 1, risk 2 = 0  

#################
dat<-data.frame(t,delw,z1,z2)
dat <- dat[which(dat$z2>0),] 
dat<-dat[order(dat$t),] # sort the data in t

## compute fitted q_j, integrated CSHs and overall survival for all (t_i,z_i) using the estimated parameters.

eLamb <- matrix(0, ncol = 2 , nrow = length(dat[,1]))

eq1 <- ea1*exp(dat$z1*eb11+dat$z2*eb12)/(ea1*exp(dat$z1*eb11+dat$z2*eb12)+ea1*ea2*exp(dat$z1*eb21+dat$z2*eb22))
eq2 <- ea1*ea2*exp(dat$z1*eb21+dat$z2*eb22)/(ea1*exp(dat$z1*eb11+dat$z2*eb12)+ea1*ea2*exp(dat$z1*eb21+dat$z2*eb22))

#semiparametric CSHs

varsigma <- matrix(0, ncol = 1 , nrow = length(dat[,1]))
for(i in 1:length(bcumh[,1])){
  I0<-which(dat$t==bcumh[i,2])
  varsigma[I0]<-bcumh[i,1]
}
eLamb[,1] <- varsigma*exp(dat$z1*eb11+dat$z2*eb12)
eLamb[,2] <- ea2*varsigma*exp(dat$z1*eb21+dat$z2*eb22)

ePi <- exp(-(eLamb[,1]+eLamb[,2]))

#estimate tau

opttau <- function(etau){
  # compute the estimated CGE
  eSCGE <- invpsigumbel(eq1*psigumbel(ePi,etau),etau)
  
 
  # estimate using regression

  y <- log(dat$t) # transform t to log(t)
  y[which(is.infinite(y))]<-NaN
  
  if (model == "expo") {
    lcge<- log(-log(eSCGE)) # transform ecge to log(-log)
    lcge[which(is.infinite(lcge))]<-NaN
    y<-y-lcge
    z1 <- dat$z1
    z2 <- dat$z2
    datlm<-data.frame(y,z1,z2) # run regression need a data frame
    fit<-lm(formula =y ~ z1+z2, datlm ) # regress log(t)-log(-log(cge)) on z, nn constant
    ea<-exp(-fit$coefficients[[1]])
    eb1<- -fit$coefficients[[2]]
    eb2<- -fit$coefficients[[3]]
    eSPAR <- exp(-dat$t*ea*exp(dat$z1*eb1+dat$z2*eb2))
  }
  if (model == "weib") {
    lcge<- log(-log(eSCGE)) # transform ecge to log(-log)
    lcge[which(is.infinite(lcge))]<-NaN
    z1 <- dat$z1
    z2 <- dat$z2
    datlm<-data.frame(y,lcge,z1,z2) # run regression need a dataframe
    fit<-lm(formula =y ~ z1+z2+lcge, datlm ) # regress log(t) on log(-log(cge)), z, ann constant
    eb <- 1/fit$coefficients[[4]] # parameter transformation
    ea<-exp(-fit$coefficients[[1]])  # parameter transformation
    eb1<- -fit$coefficients[[2]]*eb  # parameter transformation
    eb2<- -fit$coefficients[[3]]*eb  # parameter transformation
    eSPAR <- exp(-(dat$t*ea)^eb*exp(dat$z1*eb1+dat$z2*eb2))  # construct Sj implied by the param. model
  }
  
  CvM <- sum((eSCGE-eSPAR)^2) #  objective function is the mean squared sum error
  return(CvM)
}

#Set range of tau for cge for numerical minimisation of the variance criterion in theta:
c_l=0
c_u=0.95

#Set range of tau for cge for OPTIONAL plots:
taugd <-c(seq(c_l,c_u,0.05))

#determine start value for numerical minimisation

f<- rep(NA,length(taugd))
for(i in 1:length(taugd)){
  f[i] <- opttau(taugd[i])
}

I<-which(f==min(f))
startvalue<-taugd[I]

f1<-optim(startvalue, opttau, method=c("Brent"), lower=c_l,upper=c_u)
etau = f1$par

############################################
### plot of objective function
############################################


#plot objective function for assumed copula & marginal
pdf("PATH\\cvm_criterion.pdf")
plot(taugd, f, type='l',axes=FALSE, xlab=TeX(r'($\tau$)'), ylab="Criterion")
axis(1, taugd)
axis(2)
dev.off()

#estimate parameters of parametric model
#Plot estimated marginals: parametric and Cox model

#obtain Cox estimates
#risk 1
D1<-Surv(dat$t,dat$del==1)
cox_res1<-coxph(D1 ~  dat$z1 + dat$z2) # add covariates if necessary
S1_cox<-survfit(cox_res1) # estimated conditional survival, default: sample mean of covariate
S1_cox_res <- data.frame(S1_cox$time, S1_cox$surv)
S1_cox_res<-S1_cox_res[which(S1_cox_res$S1_cox.time<=400),] 


#risk 2
D2<-Surv(dat$t,dat$del==2)
cox_res2<-coxph(D2 ~  dat$z1 + dat$z2) # add covariates if necessary
S2_cox<-survfit(cox_res2) # estimated conditional survival, default: sample mean of covariate
S2_cox_res <- data.frame(S2_cox$time, S2_cox$surv)
S2_cox_res<-S2_cox_res[which(S2_cox_res$S2_cox.time<=400),]

#obtain estimated survivals of exponential model, on time grid of Cox model
mz1<-mean(dat$z1) #mean of z1
mz2<-mean(dat$z2) #mean of z2


#obtain parameters for risk 1 & 2

eS1CGE <- invpsigumbel(eq1*psigumbel(ePi,etau),etau)
eS2CGE <- invpsigumbel(eq2*psigumbel(ePi,etau),etau)
lt <- log(dat$t) # transform t to log(t)
lt[which(is.infinite(lt))]<-NaN

if (model == "expo") {
  lcge1<- log(-log(eS1CGE)) # transform ecge to log(-log)
  lcge1[which(is.infinite(lcge1))]<-NaN
  lt1<-lt-lcge1
  z1 <- dat$z1
  z2 <- dat$z2
  dat1lm<-data.frame(lt1,z1,z2) # run regression need a data frame
  fit1<-lm(formula =lt1 ~ z1+z2, dat1lm ) # regress log(t)-log(-log(cge)) on z, nn constant
  ea<- exp(-fit1$coefficients[[1]])
  eb1<- -fit1$coefficients[[2]]
  eb2<- -fit1$coefficients[[3]]
  eS1PAR <- exp(-S1_cox_res$S1_cox.time*ea*exp(mz1*eb1+mz2*eb2)) #at mean z for plot
  
  lcge2<- log(-log(eS2CGE)) # transform ecge to log(-log)
  lcge2[which(is.infinite(lcge2))]<-NaN
  lt2<-lt-lcge2
  z1 <- dat$z1
  z2 <- dat$z2
  dat2lm<-data.frame(lt2,z1,z2) # run regression need a data frame
  fit2<-lm(formula =lt2 ~ z1+z2, dat2lm ) # regress log(t)-log(-log(cge)) on z, nn constant
  ea_2<- exp(-fit2$coefficients[[1]])
  eb1_2<- -fit2$coefficients[[2]]
  eb2_2<- -fit2$coefficients[[3]]
  eS1PAR <- exp(-S2_cox_res$S2_cox.time*ea_2*exp(mz1*eb1_2+mz2*eb2_2)) #at mean z for plot
  resAFT <- list(tau=etau, a11=ea, b11=eb1, b12=eb2, a21=ea_2, b21=eb1_2, b22=eb2_2)
  #notation in paper: (tau,a_11,b_11,b_12,a_21,b_21,b_22)
}
if (model == "weib") {
  lcge1<- log(-log(eS1CGE)) # transform ecge to log(-log)
  lcge1[which(is.infinite(lcge1))]<-NaN
  z1 <- dat$z1
  z2 <- dat$z2
  dat1lm<-data.frame(lt,lcge1,z1,z2) # run regression need a dataframe
  fit1<-lm(formula =lt ~ z1+z2+lcge1, dat1lm ) # regress log(t) on log(-log(cge)), z, ann constant
  eb <- 1/fit1$coefficients[[4]] # parameter transformation
  ea<-exp(-fit1$coefficients[[1]])  # parameter transformation
  eb1<- -fit1$coefficients[[2]]*eb  # parameter transformation
  eb2<- -fit1$coefficients[[3]]*eb  # parameter transformation
  eS1PAR <- exp(-(S1_cox_res$S1_cox.time*ea)^eb*exp(mz1*eb1+mz2*eb2))  # construct Sj implied by the param. model at mean z for plot
  
  lcge2<- log(-log(eS2CGE)) # transform ecge to log(-log)
  lcge2[which(is.infinite(lcge2))]<-NaN
  z1 <- dat$z1
  z2 <- dat$z2
  dat2lm<-data.frame(lt,lcge2,z1,z2) # run regression need a dataframe
  fit2<-lm(formula =lt ~ z1+z2+lcge2, dat2lm ) # regress log(t) on log(-log(cge)), z, ann constant
  eb_2 <- 1/fit2$coefficients[[4]] # parameter transformation
  ea_2<-exp(-fit2$coefficients[[1]])  # parameter transformation
  eb1_2<- -fit2$coefficients[[2]]*eb_2  # parameter transformation
  eb2_2<- -fit2$coefficients[[3]]*eb_2  # parameter transformation
  eS2PAR <- exp(-(S2_cox_res$S2_cox.time*ea_2)^eb_2*exp(mz1*eb1_2+mz2*eb2_2))  # construct Sj implied by the param. model at mean z for plot
  resAFT <- list(tau=etau,a11=ea, a12=eb, b11=eb1, b12=eb2, a21=ea_2, a22=eb_2, b21=eb1_2, b22=eb2_2)
  #notation in paper: (tau,a_11,a_12,b_11,b_12,a_21,a_22,b_21,b_22)
}

resAFT

c(resAFT$tau, resAFT$a11, resAFT$a12, exp(resAFT$b11),exp(resAFT$b12), resAFT$a21, resAFT$a22, exp(resAFT$b21), exp(resAFT$b22))


# plot the estimated survivals
#Risk1
pdf("PATH\\application\\S1.pdf")
matplot(S1_cox_res$S1_cox.time,eS1PAR, type= "l", main="Risk 1", xlab="employment duration (in weeks)", ylab=TeX(r'($S_1(t|z)$)'), ylim=c(0, 1))
lines(S1_cox_res$S1_cox.time,S1_cox_res$S1_cox.surv, lty=2, col ="black")
legend("topright", legend =c(expression("Parametric model"(hat(tau))),expression(paste("Cox model(", tau , "=0)"))),lty=c(1,2),col=c("black","black"))                                                                                                                              
dev.off()
while (!is.null(dev.list()))  dev.off()

#Risk1
pdf("PATH\\application\\S2.pdf")
matplot(S2_cox_res$S2_cox.time,eS2PAR, type= "l", main="Risk 2", xlab="employment duration (in weeks)", ylab=TeX(r'($S_2(t|z)$)'), ylim=c(0, 1))
lines(S2_cox_res$S2_cox.time,S2_cox_res$S2_cox.surv, lty=2, col ="black")
legend("topright", legend =c(expression("Parametric model"(hat(tau))),expression(paste("Cox model(", tau , "=0)"))),lty=c(1,2),col=c("black","black"))                                                                                                                              
dev.off()
while (!is.null(dev.list()))  dev.off()


eLamb <- matrix(0, ncol = 2 , nrow = length(bcumh[,1]))
eq1 <- ea1*exp(mz1*eb11+mz2*eb12)/(ea1*exp(mz1*eb11+mz2*eb12)+ea1*ea2*exp(mz1*eb21+mz2*eb22))
eq2 <- ea1*ea2*exp(mz1*eb21+mz2*eb22)/(ea1*exp(mz1*eb11+mz2*eb12)+ea1*ea2*exp(mz1*eb21+mz2*eb22))

#compute eCGE at sample mean of z
eLamb[,1] <- bcumh[,1]*exp(mz1*eb11+mz2*eb12)
eLamb[,2] <- ea2*bcumh[,1]*exp(mz1*eb21+mz2*eb22)

ePi <- exp(-(eLamb[,1]+eLamb[,2]))

eS1CGE <- invpsigumbel(eq1*psigumbel(ePi,etau),etau)
eS2CGE <- invpsigumbel(eq2*psigumbel(ePi,etau),etau)

#Risk 1
pdf("PATH\\application\\S1_CGE.pdf")
matplot(S1_cox_res$S1_cox.time,eS1PAR, type= "l", main="Risk 1", xlab="employment duration (in weeks)", ylab=TeX(r'($S_1(t|z)$)'), ylim=c(0, 1))
lines(bcumh[,2],eS1CGE, lty=2, col ="black")
#lines(km0$time, km0$surv, lty=4, col="grey")
legend("topright", legend =c(expression("Parametric model"(hat(tau))),expression(paste("Semiparametric CGE(", hat(tau) , ")"))),lty=c(1,2),col=c("black","black"))                                                                                                                              
dev.off()
while (!is.null(dev.list()))  dev.off()

#Risk 2
pdf("PATH\\application\\S2_CGE.pdf")
matplot(S2_cox_res$S2_cox.time,eS2PAR, type= "l", main="Risk 2", xlab="employment duration (in weeks)", ylab=TeX(r'($S_2(t|z)$)'), ylim=c(0, 1))
lines(bcumh[,2],eS2CGE, lty=2, col ="black")
#lines(km0$time, km0$surv, lty=4, col="grey")
legend("topright", legend =c(expression("Parametric model"(hat(tau))),expression(paste("Semiparametric CGE(", hat(tau) , ")"))),lty=c(1,2),col=c("black","black"))                                                                                                                              
dev.off()
while (!is.null(dev.list()))  dev.off()


###########################################
# Run nonparametric bootstrap to obtain SE
###########################################

rep = 500  # number of bootstrap resamples 

bootAFT <- function(){
  sel <- 1:length(data0$wdur)
  M<- sample(sel, replace=TRUE)
  data0b <- data0[M,]
  summary(data0b)
  
  
  #################
  x<-data0b$wdur # input the name of duration variable 
  z<-cbind(rep(1,length(x)),data0b$jobno2, data0b$twdur) # input the name of covariate
  del = 2-data0b$single  # input the name of risk indicator risk 1 = 1, risk 2 = 0  
  
  #################
  y<-data.frame(x,del,z)
  
  ##################################################
  # Data reorganisation
  # replicate the data: need a y data frame y<-data.frame(x,del,z),
  # x is the duration
  # z are covariates without constant
  # risk indicator del = 1 (risk 1) or 2 (risk 2) , no censoring
  
  y1<-y[which(y[,2]==1),] 
  y2<-y[which(y[,2]==2),] 
  
  z11<-cbind(y1[,3:ncol(y1)],matrix(0,nrow=nrow(y1), ncol =ncol(y1)-2))
  z12<-cbind(matrix(0,nrow=nrow(y1), ncol =ncol(y1)-2),y1[,3:ncol(y1)])
  z21<-cbind(y2[,3:ncol(y2)],matrix(0,nrow=nrow(y2), ncol =ncol(y2)-2))
  z22<-cbind(matrix(0,nrow=nrow(y2), ncol =ncol(y2)-2),y2[,3:ncol(y2)])
  
  y11<-cbind(y1[,1],rep(1,nrow(y1)),z11)
  y12<-cbind(y1[,1],rep(0,nrow(y1)),z12)
  y21<-cbind(y2[,1],rep(0,nrow(y2)),z21)
  y22<-cbind(y2[,1],rep(1,nrow(y2)),z22)
  
  colnames(y11)<-c('X','del','z0','z1')
  colnames(y12)<-c('X','del','z0','z1')
  colnames(y21)<-c('X','del','z0','z1')
  colnames(y22)<-c('X','del','z0','z1')
  
  Y<- rbind(y11,y12,y21,y22)
  colnames(Y) <- c('x','del','a1','b11','b12', 'a2','b21','b22')
  
  #### this uses the coxph to estimate the partial likelihood estimates
  
  res<-coxph(Surv(Y[,1],Y[,2]) ~  Y$b11 + Y$b12 +Y$a2 +Y$b21+Y$b22) # add covariates if necessary
  res
  #c(b1,a2,b2) # it is to double check with the true data
  
  
  Y2<- cbind(Y[,1],Y[,2],Y[,3],Y[,4],Y[5]/52,Y[,6],Y[7],Y[8]/52)
  
  colnames(Y2) <- c('x','del','a1','b11','b12', 'a2','b21','b22')
  
  res2<-coxph(Surv(Y2[,1],Y2[,2]) ~  Y2$b11 + Y2$b12 +Y2$a2 +Y2$b21+Y2$b22) # add covariates if necessary
  res2
  
  #test<-cox.zph(res) # it is the schoenfeld test
  #print(test) # Ho: linar, p>0.05 means not rejected.
  
  #### estimate alpha 1 - the constant
  
  bcumh <- basehaz(res,centered=FALSE) # it can be use to check the baseline cum. hazard 
  fit<-lm(formula = bcumh[which(bcumh[,2]<600),1] ~ -1 + bcumh[which(bcumh[,2]<600),2] ) # regress integrated varsigma(t) on t, nn constant
  #summary(fit)
  
  bh<-fit$coefficients[[1]] # estimate the baseline cumulative hazard
  
  #plot(bcumh[which(bcumh[,2]<600),2],bcumh[which(bcumh[,2]<600),1], type="l", lty=1, ylab="cumulative hazard", xlab="Job duration") 
  #abline(a=0, b=bh, lty=2)
  #legend("bottom",legend=c("estimated baseline cumulative hazard", "fit for exponential hazard"), lty=1:2)
  
  ea1<-fit$coefficients[[1]]  #### estimated alpha 1
  ea1
  
  eb11<-res2$coefficients[[1]] # estiamted beta11
  eb12<-res2$coefficients[[2]]  # estiamted beta12
  ea2<-exp(res2$coefficients[[3]])   # estiamted alpha2
  eb21<-res2$coefficients[[4]]   # estiamted beta21
  eb22<-res2$coefficients[[5]]   # estiamted beta22
  
  
  #data preparation of original sample
  t<-data0b$wdur # input the name of duration variable 
  z1<-data0b$jobno2 #covariates
  z2<- data0b$twdur/52 
  delw = 2-data0b$single  # input the name of risk indicator risk 1 = 1, risk 2 = 0  
  
  #################
  dat<-data.frame(t,delw,z1,z2)
  dat <- dat[which(dat$z2>0),] 
  dat<-dat[order(dat$t),] # sort the data in t
  
  
  ### START OPTIM and search for tau ###
  # we need Qj and to estimate CGE.
  
  ## compute fitted q_j, integrated CSHs and overall survival for all (t_i,z_i) using the estimated parameters.
  
  eLamb <- matrix(0, ncol = 2 , nrow = length(dat[,1]))
  
  eq1 <- ea1*exp(dat$z1*eb11+dat$z2*eb12)/(ea1*exp(dat$z1*eb11+dat$z2*eb12)+ea1*ea2*exp(dat$z1*eb21+dat$z2*eb22))
  eq2 <- ea1*ea2*exp(dat$z1*eb21+dat$z2*eb22)/(ea1*exp(dat$z1*eb11+dat$z2*eb12)+ea1*ea2*exp(dat$z1*eb21+dat$z2*eb22))
  
  #semiparametric CSHs
  
  varsigma <- matrix(0, ncol = 1 , nrow = length(dat[,1]))
  for(i in 1:length(bcumh[,1])){
    I0<-which(dat$t==bcumh[i,2])
    varsigma[I0]<-bcumh[i,1]
  }
  eLamb[,1] <- varsigma*exp(dat$z1*eb11+dat$z2*eb12)
  eLamb[,2] <- ea2*varsigma*exp(dat$z1*eb21+dat$z2*eb22)
  
  ePi <- exp(-(eLamb[,1]+eLamb[,2]))
  
  #eSCGE <- invpsigumbel(eq1*psigumbel(ePi,0),0)
  
  #estimate tau
  
  opttau <- function(etau){
    # compute the estimated CGE
    eSCGE <- invpsigumbel(eq1*psigumbel(ePi,etau),etau)
    
    
    # estimate alpha using regression
    
    y <- log(dat$t) # transform t to log(t)
    y[which(is.infinite(y))]<-NaN
    
    if (model == "expo") {
      lcge<- log(-log(eSCGE)) # transform ecge to log(-log)
      lcge[which(is.infinite(lcge))]<-NaN
      y<-y-lcge
      z1 <- dat$z1
      z2 <- dat$z2
      datlm<-data.frame(y,z1,z2) # run regression need a data frame
      fit<-lm(formula =y ~ z1+z2, datlm ) # regress log(t)-log(-log(cge)) on z, nn constant
      ea<-exp(-fit$coefficients[[1]])
      eb1<- -fit$coefficients[[2]]
      eb2<- -fit$coefficients[[3]]
      eSPAR <- exp(-dat$t*ea*exp(dat$z1*eb1+dat$z2*eb2))
    }
    if (model == "weib") {
      lcge<- log(-log(eSCGE)) # transform ecge to log(-log)
      lcge[which(is.infinite(lcge))]<-NaN
      z1 <- dat$z1
      z2 <- dat$z2
      datlm<-data.frame(y,lcge,z1,z2) # run regression need a dataframe
      fit<-lm(formula =y ~ z1+z2+lcge, datlm ) # regress log(t) on log(-log(cge)), z, ann constant
      eb <- 1/fit$coefficients[[4]] # parameter transformation
      ea<-exp(-fit$coefficients[[1]])  # parameter transformation
      eb1<- -fit$coefficients[[2]]*eb  # parameter transformation
      eb2<- -fit$coefficients[[3]]*eb  # parameter transformation
      eSPAR <- exp(-(dat$t*ea)^eb*exp(dat$z1*eb1+dat$z2*eb2))  # construct Sj implied by the param. model
    }
    
    CvM <- sum((eSCGE-eSPAR)^2) #  objective function is the mean squared sum error
    return(CvM)
  }
  
  #Set range of tau for cge for numerical minimisation of the variance criterion in theta:
  c_l=0
  c_u=0.95
  
  #Set range of tau for cge for OPTIONAL plots:
  taugd <-c(seq(c_l,c_u,0.05))
  
  #determine start value for numerical minimisation
  
  f<- rep(NA,length(taugd))
  for(i in 1:length(taugd)){
    f[i] <- opttau(taugd[i])
  }
  
  I<-which(f==min(f))
  startvalue<-taugd[I]
  
  f1<-optim(startvalue, opttau, method=c("Brent"), lower=c_l,upper=c_u)
  etau = f1$par
  
  #estimate parameters of parametric model
  #Plot estimated marginals: parametric and Cox model
  
  #obtain Cox estimates
  #risk 1
  D1<-Surv(dat$t,dat$del==1)
  cox_res1<-coxph(D1 ~  dat$z1 + dat$z2) # add covariates if necessary
  S1_cox<-survfit(cox_res1) # estimated conditional survival, default: sample mean of covariate
  S1_cox_res <- data.frame(S1_cox$time, S1_cox$surv)
  S1_cox_res<-S1_cox_res[which(S1_cox_res$S1_cox.time<=400),] 
  
  
  #risk 2
  D2<-Surv(dat$t,dat$del==2)
  cox_res2<-coxph(D2 ~  dat$z1 + dat$z2) # add covariates if necessary
  S2_cox<-survfit(cox_res2) # estimated conditional survival, default: sample mean of covariate
  S2_cox_res <- data.frame(S2_cox$time, S2_cox$surv)
  S2_cox_res<-S2_cox_res[which(S2_cox_res$S2_cox.time<=400),]
  
  #obtain estimated survivals of exponential model, on time grid of Cox model
  mz1<-mean(dat$z1) #mean of z1
  mz2<-mean(dat$z2) #mean of z2
  
  
  #obtain parameters for risk 1 & 2
  
  eS1CGE <- invpsigumbel(eq1*psigumbel(ePi,etau),etau)
  eS2CGE <- invpsigumbel(eq2*psigumbel(ePi,etau),etau)
  lt <- log(dat$t) # transform t to log(t)
  lt[which(is.infinite(lt))]<-NaN
  
  if (model == "expo") {
    lcge1<- log(-log(eS1CGE)) # transform ecge to log(-log)
    lcge1[which(is.infinite(lcge1))]<-NaN
    lt1<-lt-lcge1
    z1 <- dat$z1
    z2 <- dat$z2
    dat1lm<-data.frame(lt1,z1,z2) # run regression need a data frame
    fit1<-lm(formula =lt1 ~ z1+z2, dat1lm ) # regress log(t)-log(-log(cge)) on z, nn constant
    ea<- exp(-fit1$coefficients[[1]])
    eb1<- -fit1$coefficients[[2]]
    eb2<- -fit1$coefficients[[3]]
    eS1PAR <- exp(-S1_cox_res$S1_cox.time*ea*exp(mz1*eb1+mz2*eb2)) #at mean z for plot
    
    lcge2<- log(-log(eS2CGE)) # transform ecge to log(-log)
    lcge2[which(is.infinite(lcge2))]<-NaN
    lt2<-lt-lcge2
    z1 <- dat$z1
    z2 <- dat$z2
    dat2lm<-data.frame(lt2,z1,z2) # run regression need a data frame
    fit2<-lm(formula =lt2 ~ z1+z2, dat2lm ) # regress log(t)-log(-log(cge)) on z, nn constant
    ea_2<- exp(-fit2$coefficients[[1]])
    eb1_2<- -fit2$coefficients[[2]]
    eb2_2<- -fit2$coefficients[[3]]
    eS1PAR <- exp(-S2_cox_res$S2_cox.time*ea_2*exp(mz1*eb1_2+mz2*eb2_2)) #at mean z for plot
    #list(etau, ea, eb1, eb2, ea_2, eb1_2, eb2_2)
    #notation in paper: (tau,a_11,b_11,b_12,a_21,b_21,b_22)
  }
  if (model == "weib") {
    lcge1<- log(-log(eS1CGE)) # transform ecge to log(-log)
    lcge1[which(is.infinite(lcge1))]<-NaN
    z1 <- dat$z1
    z2 <- dat$z2
    dat1lm<-data.frame(lt,lcge1,z1,z2) # run regression need a dataframe
    fit1<-lm(formula =lt ~ z1+z2+lcge1, dat1lm ) # regress log(t) on log(-log(cge)), z, ann constant
    eb <- 1/fit1$coefficients[[4]] # parameter transformation
    ea<-exp(-fit1$coefficients[[1]])  # parameter transformation
    eb1<- -fit1$coefficients[[2]]*eb  # parameter transformation
    eb2<- -fit1$coefficients[[3]]*eb  # parameter transformation
    eS1PAR <- exp(-(S1_cox_res$S1_cox.time*ea)^eb*exp(mz1*eb1+mz2*eb2))  # construct Sj implied by the param. model at mean z for plot
    
    lcge2<- log(-log(eS2CGE)) # transform ecge to log(-log)
    lcge2[which(is.infinite(lcge2))]<-NaN
    z1 <- dat$z1
    z2 <- dat$z2
    dat2lm<-data.frame(lt,lcge2,z1,z2) # run regression need a dataframe
    fit2<-lm(formula =lt ~ z1+z2+lcge2, dat2lm ) # regress log(t) on log(-log(cge)), z, ann constant
    eb_2 <- 1/fit2$coefficients[[4]] # parameter transformation
    ea_2<-exp(-fit2$coefficients[[1]])  # parameter transformation
    eb1_2<- -fit2$coefficients[[2]]*eb_2  # parameter transformation
    eb2_2<- -fit2$coefficients[[3]]*eb_2  # parameter transformation
    eS2PAR <- exp(-(S2_cox_res$S2_cox.time*ea_2)^eb_2*exp(mz1*eb1_2+mz2*eb2_2))  # construct Sj implied by the param. model at mean z for plot
    #list(etau,ea, eb, eb1, eb2, ea_2, eb_2, eb1_2, eb2_2)
    #notation in paper: (tau,a_11,a_12,b_11,b_12,a_21,a_22,b_21,b_22)
  }
  
  
  if (model != "expo") { 
    res <- list(etaub = etau, ea11b= ea, ea12b = eb, eb11b=eb1, eb12b=eb2, ea21b = ea_2, ea22b = eb_2, eb21b=eb1_2, eb22b=eb2_2)
  }
  if (model == "expo") { 
    res <- list(etaub = etau, ea11b= ea, eb11b=eb1, eb12b=eb2, ea21b = ea_2, eb21b=eb1_2, eb22b=eb2_2)
  }
  
  return(res)
}


res2 <- vector("list", rep)
for(i in 1:rep) res2[[i]] <- try(bootAFT(), TRUE)
res3<-res2[sapply(res2, function(x) !inherits(x, "try-error"))]

etaub<-rep(0,length(res3)) # create data frame for estimated tau
ea11b<-rep(0,length(res3)) # create data frame for estimated beta
ea12b<-rep(0,length(res3)) # create data frame for estimated beta
eb11b<-rep(0,length(res3)) # create data frame for estimated beta
eb12b<-rep(0,length(res3)) # create data frame for estimated beta
ea21b<-rep(0,length(res3)) # create data frame for estimated beta
ea22b<-rep(0,length(res3)) # create data frame for estimated beta
if (model != "expo") { 
eb21b<-rep(0,length(res3)) # create data frame for estimated beta
eb22b<-rep(0,length(res3)) # create data frame for estimated beta
}

if (model!="expo"){
  for (i in 1:length(res3)){
    etaub[i]<-res3[[i]][[1]]
    ea11b[i]<-res3[[i]][[2]]  
    ea12b[i]<-res3[[i]][[3]]
    eb11b[i]<-res3[[i]][[4]]
    eb12b[i]<-res3[[i]][[5]]
    ea21b[i]<-res3[[i]][[6]]  
    ea22b[i]<-res3[[i]][[7]]
    eb21b[i]<-res3[[i]][[8]]
    eb22b[i]<-res3[[i]][[9]]
  }
}
if (model=="expo"){
  for (i in 1:length(res3)){
    etaub[i]<-res3[[i]][[1]]
    ea11b[i]<-res3[[i]][[2]]  
    eb11b[i]<-res3[[i]][[3]]
    eb12b[i]<-res3[[i]][[4]]
    ea21b[i]<-res3[[i]][[5]]  
    eb21b[i]<-res3[[i]][[6]]
    eb22b[i]<-res3[[i]][[7]]
  }
} 

mtau<-mean(etaub)
ma11<-mean(ea11b)
mb11<-mean(eb11b)
mb12<-mean(eb12b)
ma21<-mean(ea21b)
mb21<-mean(eb21b)
mb22<-mean(eb22b)
if (model!="expo"){
ma12<-mean(ea12b)
ma22<-mean(ea22b)
}

setau <- sqrt(mean((etaub-mtau)^2))
sea11 <- sqrt(mean((ea11b-ma11)^2))
seb11 <- sqrt(mean((eb11b-mb11)^2))
seb12 <- sqrt(mean((eb12b-mb12)^2))
sea21 <- sqrt(mean((ea21b-ma21)^2))
seb21 <- sqrt(mean((eb21b-mb21)^2))
seb22 <- sqrt(mean((eb22b-mb22)^2))
if (model!="expo"){
sea12 <- sqrt(mean((ea12b-ma12)^2))
sea22 <- sqrt(mean((ea22b-ma22)^2))
}

#standard errors
setau # SE for tau
sea11 # SE for a11
seb11 # SE for b11
seb12 # SE for b12
sea21 # SE for a21
seb21 # SE for b21
seb22 # SE for b22
if (model!="expo"){
  sea12 # SE for a12
}
if (model!="expo"){
  sea22 # SE for a22
}

#p-values
z_tau<-(resAFT$tau)/setau
2*pnorm(abs(z_tau), lower.tail=FALSE)
z_a11<-(resAFT$a11)/sea11
2*pnorm(abs(z_a11), lower.tail=FALSE)
z_b11<-(resAFT$b11)/seb11
2*pnorm(abs(z_b11), lower.tail=FALSE)
z_b12<-(resAFT$b12)/seb12
2*pnorm(abs(z_b12), lower.tail=FALSE)
z_a21<-(resAFT$a21)/sea21
2*pnorm(abs(z_a21), lower.tail=FALSE)
z_b21<-(resAFT$b21)/seb21
2*pnorm(abs(z_b21), lower.tail=FALSE)
z_b22<-(resAFT$b22)/seb22
2*pnorm(abs(z_b22), lower.tail=FALSE)
if (model!="expo"){
  z_a12<-(resAFT$a12-1)/sea12
  2*pnorm(abs(z_a12), lower.tail=FALSE)
}
if (model!="expo"){
  z_a22<-(resAFT$a22-1)/sea22
  2*pnorm(abs(z_a22), lower.tail=FALSE)
}