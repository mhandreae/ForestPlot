Ware = read.csv('WarePostTreat.csv',header=T,na.str=c(".","@@"))
outcome = sapply(split(Ware[Ware[,4]<6,7],list(Ware[Ware[,4]<6,1],Ware[Ware[,4]<6,3])),mean,na.rm=T)
conditions = strsplit(names(outcome),".",fixed=T)
dat = cbind(outcome, matrix(as.integer(unlist(conditions)),length(outcome),2,byrow=T))
treatord = read.csv('WareTreatOrder.csv',header=T)[1:23,1:2]
trt = strsplit(as.character(treatord[,1]),"")
names(trt) = treatord[,2]
trt = as.data.frame(t(as.data.frame(trt)))
subj = strsplit(rownames(trt),"X")
subj = (lapply(subj, function(v){v[2]}))
subj = unlist(subj)
rownames(trt) = subj
tt = as.numeric(apply(dat[,2:3],1,function(v){trt[as.character(v[1]),v[2]]}))
dat = cbind(dat,tt)
dat[,2] = as.numeric(as.factor(dat[,2]))
base = read.csv('WareBase.csv',header=T)
base[,1] = as.numeric(as.factor(base[,1]))
dat = cbind(dat,baseline=base[dat[,2],2])
colnames(dat)= c('DayAvg','ID','Cycle','Treat','baseline')
dat = cbind(dat, resp=ifelse(dat[,1]<.7*dat[,5],1,0))
mod1 = glmer(resp ~ factor(Treat) + (1|ID) , family=binomial, link='logit',data=as.data.frame(dat))
mod0 = glmer(resp ~ (1|ID) , family = binomial, link='logit',data=as.data.frame(dat))
anova(mod0,mod1)
dose = c(0,2.5,6,9.4)
dat = cbind(dat, dose = dose[dat[,4]])
Ware.dat = dat

################

Ellis = matrix(c(1,1,1,0,0,1,0,0),4,2,byrow=T)[c(rep(1,3),rep(2,2),rep(3,10),rep(4,13)),]
colnames(Ellis) = c('d0','d96')
Ellis.dat = cbind(ID=rep(1:28,2),Treat=rep(1:2,each=28),resp=as.vector(Ellis),dose=rep(c(0,96),each=28))

#################

Wilsey = matrix(c(
0,0,0,
1,1,1,
1,1,1,
0,0,0,
1,1,1,
0,1,0,
1,1,0,
1,1,1,
0,0,1,
1,1,1,
1,NA,NA,
1,1,1,
0,NA,NA,
1,1,0,
0,1,0,
0,0,0,
1,1,0,
1,1,1,
1,1,1,
1,1,1,
0,0,0,
0,0,1,
0,1,0,
1,1,1,
1,0,0,
1,0,0,
1,1,1,
1,1,0,
1,1,1,
1,1,0,
1,1,1,
NA,1,1,
0,0,0,
1,1,1,
NA,0,NA,
NA,0,1,
NA,0,NA,
NA,1,NA),38,3, byrow=T)
colnames(Wilsey) = c('d34.3','d19.25','d0')
Wilsey.dat = data.frame(ID = rep(1:38,3), Treat = rep(3:1,each=38), resp = as.vector(Wilsey), dose = rep(c(34.3,19.25,0),each=38))


#################

Wilsey2 =matrix(c(1,1,1,
0,0,0,
0,0,0,
0,0,0,
0,0,0,
1,0,0,
1,0,0,
1,0,1,
0,0,0,
1,1,0,
0,0,0,
0,0,0,
1,1,0,
0,0,0,
1,1,1,
0,0,0,
0,0,0,
1,0,0,
0,0,0,
1,1,0,
1,1,0,
1,1,1,
0,0,0,
1,1,0,
1,1,0,
0,1,0,
1,0,1,
0,1,1,
1,0,1,
1,1,0,
0,0,0,
1,1,1,
0,0,1,
0,NA,NA,
0,1,1,
1,1,0,
NA,1,0,
NA,1,1,
NA,NA,0), 39,3,byrow=T)
colnames(Wilsey2) = c('d18','d9','d0')
Wilsey2.dat = data.frame(ID = rep(1:39,3), Treat = rep(3:1,each=39), resp = as.vector(Wilsey2), dose = rep(c(18,9,0), each=39))

#####################

Abrams.dat = data.frame(ID=1:50, Treat = rep(1:2, each=25), resp = c(rep(1,6),rep(0,19),rep(1,13),rep(0,12)), dose = rep(c(0,96),each=25))

########################

AllResp = rbind(Abrams.dat, Ellis.dat, Ware.dat[,c(2,4,6,7)], Wilsey.dat,Wilsey2.dat)
AllResp = cbind(AllResp, Study = c(rep(1,nrow(Abrams.dat)), rep(2, nrow(Ellis.dat)), rep(3,nrow(Ware.dat)), rep(4,nrow(Wilsey.dat)), rep(5,nrow(Wilsey2.dat))))
AllResp$ID = factor(AllResp$ID)
AllResp$Study = factor(AllResp$Study)
AllResp$Treat = factor(AllResp$Treat)

mod.ldose = glmer(resp ~ (1|ID:Study) + log(1+dose) + (log(1+dose) | Study),data=AllResp,family=binomial(logit))
mod.dose = glmer(resp ~ (1|ID:Study) + dose + (dose | Study),data=AllResp,family=binomial(logit))
mod.ldose2 = glmer(resp ~ (1|ID:Study) + log(1+dose) + (log(1+dose) | Study) + (1|Treat:Study),data=AllResp,family=binomial(logit))
mod.ldose3 = glmer(resp ~ (1|ID:Study) + log(1+dose)  + (1|Treat:Study),data=AllResp,family=binomial(logit))

mycoef = fixef(mod.ldose) + t(ranef(mod.ldose,which='Study')$Study)
dd = seq(0,100,len=200)
x = cbind(1,log(1+dd))
psi = x%*%mycoef
matplot(dd,1/(1+exp(-psi)),type='l',col=gray(.5),lty=1,xlab='Dose',ylab='Probabilty of Response',ylim=c(0,1))
lines(dd,1/(1+exp(-x%*%fixef(mod.ldose))),lwd=2)

# individual estimates #
slopeinf = sapply(1:5, function(i){ mod = glmer(resp~log(1+dose) + (1|ID),family=binomial,data=AllResp[AllResp$Study==i,]); b = fixef(mod)[2]; v=vcov(mod)[2,2]; c(b,sqrt(v))})

metaplot(slopeinf[1,],slopeinf[2,],summn=.43185,sumse=.08596,sumnn=1/(.08596^2),labels=c('Abrams','Ellis','Ware','Wilsey 1','Wilsey 2'),xlab='Dose Effect (Slope)',main='Frequentist',ylab='Study')
dev.print(pdf,file='FreqForest.pdf',width=6,height=6)
forestplot(cbind(c('Abrams','Ellis','Ware','Wilsey 1','Wilsey 2','Frequentist','Bayesian')), c(slopeinf[1,],.43185,.4632), c(slopeinf[1,]-2*slopeinf[2,], .43185-2*.08596,.1953),c(slopeinf[1,]+2*slopeinf[2,], .43185+2*.08596,.7535),is.summary=c(F,F,F,F,F,T,T))
forestplot(cbind(c('Abrams','Ellis','Ware','Wilsey 1','Wilsey 2','','Frequentist','','Bayesian')), c(slopeinf[1,],NA,.43185,NA,.4632), c(slopeinf[1,]-2*slopeinf[2,], NA,.43185-2*.08596,NA,.1953),c(slopeinf[1,]+2*slopeinf[2,],NA, .43185+2*.08596,NA,.7535),is.summary=c(F,F,F,F,F,F,T,F,T))
dev.print(pdf,file='BothForest.pdf',width=6,height=6)
                                        # Tried the following model because correlation between slope and intercept was -1 #
mod.ldose.nocor = glmer(resp ~ (1|ID:Study) + log(1+dose) + (0+log(1+dose) | Study) + (1|Study),data=AllResp,family=binomial(logit))
mod.ldose.nocor2 = glmer(resp ~ (1|ID:Study) + log(1+dose)  + (1|Study),data=AllResp,family=binomial(logit))


forestplot(cbind(c('Abrams','Ellis','Ware','Wilsey 1','Wilsey 2','','Frequentist','Frequentist NA->0','','Bayesian')), c(1.519,1.54,1.308, 1.453,1.448,NA,1.3302,1.2335,NA,1.46), c(.5111,.5818,.08471,.5711,.5596,NA,1.3302-2*.2803,1.2335-2*.27,NA,.5945),c(2.896,2.75,2.283,2.474,2.398,NA,1.3302+2*.2803,1.2335+2*.27,NA,2.437), is.summary=c(F,F,F,F,F,F,T,T,F,T),main='Ignoring Dose',xlab='Log Odds Ratio')

forestplot(cbind(c('Study','Abrams','Ellis','Ware','Ware','Ware','Wilsey 1','Wilsey 1','Wilsey 2','Wilsey 2','','Bayesian'),c('Dose',96,96,2.5,6.3,9.4,19.25,34.3,9,18,'','')), c(NA,1.557,1.669,1.194, 1.239,1.356,1.337,1.331,1.362,1.456,NA,1.43), c(NA,.4842,.6614,-.0921,-0.0022,.1804,.3736,.3582,.4569,.5773,NA,.5863),c(NA,2.94,2.886,2.326,2.388,2.587,2.327,2.328,2.321,2.462,NA,2.401), is.summary=c(T,F,F,F,F,F,F,F,F,F,F,T),main='Ignoring Dose',xlab='Log Odds Ratio')

dev.print(pdf,file='BayesianForestNoDose.pdf',width=6,height=6)


########

summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==1 & (AllResp$Treat== 1 | AllResp$Treat==2),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==2 & (AllResp$Treat== 1 | AllResp$Treat==2),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==3 & (AllResp$Treat== 1 | AllResp$Treat==2),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==3 & (AllResp$Treat== 1 | AllResp$Treat==3),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==3 & (AllResp$Treat== 1 | AllResp$Treat==4),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==4 & (AllResp$Treat== 1 | AllResp$Treat==2),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==4 & (AllResp$Treat== 1 | AllResp$Treat==3),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==5 & (AllResp$Treat== 1 | AllResp$Treat==2),] ,family=binomial(link='logit')))
summary(glmer(resp ~ Treat + (1|ID),data=AllResp[AllResp$Study==5 & (AllResp$Treat== 1 | AllResp$Treat==3),] ,family=binomial(link='logit')))

tmp.fun = function(dat){
  mod = glmer(resp~ Treat + (1|ID), family=binomial,data=dat)
  arms = fixef(mod)[-1]
  se.arms = sqrt(diag(vcov(mod)))[-1]
  K = length(arms)
  est = mean(arms)
  se.est = sqrt(sum(vcov(mod)[-1,-1]))/K
  cbind(c(arms,est),c(se.arms,se.est))
}

myOR.fun = function(tbl){
  N = sum(tbl)
  phat = as.vector(tbl/N)
  v.phat = (diag(phat) - tcrossprod(phat))/N
  lphat = log(phat)
  v.lphat = sweep(v.phat/phat,2,phat,'/')
  est = lphat[3]-lphat[2]
  cc = c(0,1,-1,0)
  v.est = crossprod(cc,v.lphat)%*%cc
  se.est = sqrt(v.est)
  c(est,est+1.96*c(-1,1)*se.est)
}
  
  

forestplot(cbind(c('Study','Abrams','Ellis','Ware','Ware','Ware','Wilsey 1','Wilsey 1','Wilsey 2','Wilsey 2','','Bayesian',NA,NA,NA),c('Dose',96,96,2.5,6.3,9.4,19.25,34.3,9,18,'','',"","","")), c(NA,1.2327,1.609,.405, 1.099,1.609,.981,1.253,.916,1.299,NA,1.43,1.139,.765), c(NA,0,.091,-1.38,-1.165,-.538,-.346,-.319,-.243,.023,NA,.5863,.448,.332),c(NA,2.4649,3.13,2.195,3.36,3.756,2.308,2.824,2.076,2.576,NA,2.401,1.951,1.237), is.summary=c(T,F,F,F,F,F,F,F,F,F,F,T,T,T),main='Ignoring Dose',xlab='Log Odds Ratio')

out = read.table('NoDoseOut.txt')[,-1]
tmp.out = matrix(out[c(160001:180000,200001:220000,230001:310000)],20000,6)
tmp.out = read.table('tmpOut.txt')[,-1]
tmp.out = matrix(tmp.out,100000,6)
p1.fun = function(v){integrate(function(z){ dnorm(z,0,sqrt(sum(v[3:4])))/(1+exp(-v[1]-v[2]-z))}, -Inf,Inf)$val}
p0.fun = function(v){integrate(function(z){ dnorm(z,0,sqrt(v[3]))/(1+exp(-v[2]-z))}, -Inf,Inf)$val}
quantile(apply(tmp.out[1:10000,],1,function(v){p1 = p1.fun(v); p0 = p0.fun(v); log(p1*(1-p0))-log(p0*(1-p1))}),c(.5,.025,.975))
p1.fun = function(v){integrate(function(z){ dnorm(z,0,sqrt(sum(v[3:5])))/(1+exp(-v[1]-v[2]-z))}, -Inf,Inf)$val}
p0.fun = function(v){integrate(function(z){ dnorm(z,0,sqrt(v[3]+v[5]))/(1+exp(-v[2]-z))}, -Inf,Inf)$val}
quantile(apply(tmp.out[1:10000,],1,function(v){p1 = p1.fun(v); p0 = p0.fun(v); log(p1*(1-p0))-log(p0*(1-p1))}),c(.5,.025,.975))
quantile(apply(tmp.out[1:10000,],1,function(v){p1 = p1.fun(v); p0 = p0.fun(v); 1/(p1-p0)}),c(.5,.025,.975))

forestplot(cbind(c('Study','Abrams','Ellis','Ware','Ware','Ware','Wilsey 1','Wilsey 1','Wilsey 2','Wilsey 2','','Bayesian',NA,NA,NA),c('Dose',96,96,2.5,6.3,9.4,19,34,9.0,18,'','',"","",""),c("Placebo","6/25","5/28","3/22","","","18/33","","11/38","","","","","",""),c("Treat","13/25","13/28","4/21","5/22","7/21","24/36","22/33","17/37","18/36","","","","","")), c(NA,1.2327,1.609,.405, 1.099,1.609,.981,1.253,.916,1.299,NA,1.43,1.139,.765), c(NA,0,.091,-1.38,-1.165,-.538,-.346,-.319,-.243,.023,NA,.5863,.448,.332),c(NA,2.4649,3.13,2.195,3.36,3.756,2.308,2.824,2.076,2.576,NA,2.401,1.951,1.237), is.summary=c(T,F,F,F,F,F,F,F,F,F,F,T,T,T),main='Ignoring Dose',xlab='Log Odds Ratio')

gob = apply(tmp.out[1:10000,],1,function(v){p1 = p1.fun(v); p0 = p0.fun(v); log(p1*(1-p0))-log(p0*(1-p1))})
gob.wt = dnorm(tmp.out[1:10000,1],0,1e4)/dnorm(tmp.out[1:10000,1],0,10)
gob2 = cbind(gob,gob.wt/sum(gob.wt))[order(gob),]
gob2[cumsum(gob2[,2])>.975,][1:5,]
# changing prior on Beta moved CI from (.448,1.951) (median=1.14) to (.449,1.954) (median=1.141) #
gob = apply(tmp.out[1:10000,],1,function(v){p1 = p1.fun(v); p0 = p0.fun(v); log(p1*(1-p0))-log(p0*(1-p1))})
gob.wt = dt(tmp.out[1:10000,4],1)/dt(tmp.out[1:10000,4]/100,1)
gob2 = cbind(gob,gob.wt/sum(gob.wt))[order(gob),]
gob2[cumsum(gob2[,2])>.975,][1:5,]
gob2[cumsum(gob2[,2])>.025,][1:5,]
gob2[cumsum(gob2[,2])>.5,][1:5,]

sens.fun = function(sf){
  gob.wt = dt(tmp.out[1:10000,4]/sf,1)/dt(tmp.out[1:10000,4],1)
  gob2 = cbind(gob,gob.wt/sum(gob.wt))[order(gob),]
  q2 = gob2[cumsum(gob2[,2])>.975,][1,1]
  q1 = gob2[cumsum(gob2[,2])>.025,][1,1]
  m = gob2[cumsum(gob2[,2])>.5,][1,1]
  c(m,q1,q2)
}
sens.out = sapply(seq(.1,20,len=100),sens.fun)
matplot(seq(.1,20,len=100),matrix(sens.out,100,3,byrow=T),type='l',lty=c(1,2,2),col=1,xlab='Scale Factor',ylab='Log Odds Ratio')
abline(v=1)

forestplot(cbind(c('Study','','Abrams 07','','Ellis 09','','Ware 10','Ware 10','Ware 10','','Wilsey 08','Wilsey 08','','Wilsey 13','Wilsey 13','','Bayesian'),c('Dose','',96,'',96,'',2.5,6.3,9.4,'',19,34,'',9.0,18,'',''),c("Placebo",'',"6/25",'',"5/28",'',"3/22","","",'',"18/33","",'',"11/38","","",""),c("Treat",'',"13/25",'',"13/28",'',"4/21","5/22","7/21",'',"24/36","22/33",'',"17/37","18/36","",""),
           c('Est. OR (CI)',NA, '3.43 (1.00,11.8)', NA, '5.00 (1.10,22.9)',NA, '1.50 (0.25,8.98)',
             '3.00 (0.31,28.8)','5.00 (0.58,42.8)',NA,'2.67 (0.71,10.1)','3.50 (0.73,16.8)',NA,
             '2.50 (0.78,7.97)','3.67 (1.02,13.1)', NA, '3.22 (1.59,7.24)')),
           (c(NA,NA,1.2327,NA,1.609,NA,.405, 1.099,1.609,NA,.981,1.253,NA,.916,1.299,NA,1.169)), (c(NA,NA,0,NA,.091,NA,-1.38,-1.165,-.538,NA,-.346,-.319,NA,-.243,.023,NA,.462)),(c(NA,NA,2.4649,NA,3.13,NA,2.195,3.36,3.756,NA,2.308,2.824,NA,2.076,2.576,NA,1.986)), is.summary=c(T,T,F,T,F,T,F,F,F,T,F,F,T,F,F,T,T),main='Ignoring Dose',xlab='Odds Ratio',xlog=T,xlim=c(-1,exp(4)),xticks=c(.2,.5,1,2,5,10,20,40))

or.stats = round(cbind(exp(c(NA,NA,1.2327,NA,1.609,NA,.405, 1.099,1.609,NA,.981,1.253,NA,.916,1.299,NA,1.43,1.139,.765,1.168)), exp(c(NA,NA,0,NA,.091,NA,-1.38,-1.165,-.538,NA,-.346,-.319,NA,-.243,.023,NA,.5863,.448,.332,.288)),exp(c(NA,NA,2.4649,NA,3.13,NA,2.195,3.36,3.756,NA,2.308,2.824,NA,2.076,2.576,NA,2.401,1.951,1.237,2.179))),2)

