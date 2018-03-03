#Example codes to implement riBART

#Author: Yaoyuan Vincent, Tan
#Date created: 22 Feb 2018
#Date updated: 22 Feb 2018

rm(list=ls())

#####Continuous outcomes######
#Setup toy dataset
nk=5
K=50
num=nk*K
true_tau=1
sd=1
p=10
X=matrix(runif(nk*K*p),ncol=p)
G=10*sin(pi*X[,1]*X[,2])+20*(X[,3]-.5)^2+10*X[,4]+5*X[,5]
E=rnorm(num,0,sd)
A=rnorm(K,0,true_tau)
y=G+rep(A,each=nk)+E
groupings=rep(1:K,each=nk)
gindex=sort(unique(groupings))

#Setup ribart initial parameters
library(pscl)
library(zipfR)
library(mda)
numpost=5000
numburnin=1000
R_init=mars(X,y)$residuals
initak=numeric(K)
for(i in 1:K) initak[i]=mean(R_init[groupings==gindex[i]])
nu_sig=3
inv_rgamma=Rgamma.inv(nu_sig/2,0.1,lower=F)
zsig2_ribart=sum((R_init-c(t(matrix(rep(initak,nk),K,nk))))^2)/(nk*K-(nk*K*(1-sqrt(sum(R_init^2)/(mars(X,y)$gcv*nk*K)))+K))
lambda_ribart=2*zsig2_ribart*inv_rgamma/nu_sig
nu=2
lambda=1
alpha=0.95
beta=2
k=2
m=200
base=1.96
sigma_mu=base/(k*sqrt(m))
probgrow=0.5
probprune=0.5

#riBART
source("bart_fns.R")

#initializations
initsigma=sqrt(rigamma(1,(num+nu_sig)/2,(nu_sig*lambda_ribart+sum((R_init-rep(initak,each=nk))^2))/2))
inittau_prop=sqrt(rigamma(1,(K+nu)/2,(nu*lambda+sum(initak^2))/2))

initak_prop=rnorm(K)

for (i in 1:K){
	sum_R_init=sum(R_init[groupings==gindex[i]])
	initak_prop[i]=initak_prop[i]*(initsigma*inittau_prop)/sqrt(nk*inittau_prop^2+initsigma^2)+(inittau_prop^2*sum_R_init)/(nk*inittau_prop^2+initsigma^2)
}

tildey_prop=y-rep(initak_prop,each=nk)

beta1=(2*base)/(max(tildey_prop)-min(tildey_prop))
beta0=base-max(tildey_prop)*beta1
transy=beta0+beta1*tildey_prop
sigma=beta1*initsigma
init_mu_ij=matrix(1/m*mean(transy),length(transy),m)
init_value=list()
init_xindex=list()
init_depth=list()
init_varindex=list()
init_variable=list()
for(i in 1:m){
	init_value[[i]]=NA
	init_xindex[[i]]=NA
	init_depth[[i]]=1
	init_varindex[[i]]=1
	init_variable[[i]]=rep(1,length(transy))	
}
onedraw_prop=bart1draw(transy,init_mu_ij,sigma,sigma_mu,X,init_value,init_xindex,init_depth,init_varindex,init_variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
gxprop=onedraw_prop$un_trans
R_prop=y-gxprop

#burnin
for(burn in 1:numburnin){
	bsig_prop=sqrt(rigamma(1,(num+nu_sig)/2,(nu_sig*lambda_ribart+sum((tildey_prop-gxprop)^2))/2))
	inittau_prop=sqrt(rigamma(1,(K+nu)/2,(nu*lambda+sum(initak_prop^2))/2))
	initak_prop=rnorm(K)
	for (i in 1:K){
		sum_R_prop=sum(R_prop[groupings==gindex[i]])
		initak_prop[i]=initak_prop[i]*(bsig_prop*inittau_prop)/sqrt(nk*inittau_prop^2+bsig_prop^2)+(inittau_prop^2*sum_R_prop)/(nk*inittau_prop^2+bsig_prop^2)
	}
	tildey_prop=y-rep(initak_prop,each=nk) 
	beta1=(2*base)/(max(tildey_prop)-min(tildey_prop))
	beta0=base-max(tildey_prop)*beta1
	transy=beta0+beta1*tildey_prop
	sigma=beta1*bsig_prop
	onedraw_prop=bart1draw(transy,onedraw_prop$mu_ij,sigma,sigma_mu,X,onedraw_prop$value,onedraw_prop$xindex,onedraw_prop$depth,onedraw_prop$varindex,onedraw_prop$variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
	gxprop=onedraw_prop$un_trans
	R_prop=y-gxprop
}

#posterior draws of sigma, tau, gx, and the random intercept
psig_prop=numeric(numpost)
ptau_prop=numeric(numpost)
pgxprop=matrix(0,nrow=num,ncol=numpost)
pakprop=matrix(0,nrow=K,ncol=numpost)

for(post in 1:numpost){
	bsig_prop=sqrt(rigamma(1,(num+nu_sig)/2,(nu_sig*lambda_ribart+sum((tildey_prop-gxprop)^2))/2))
	inittau_prop=sqrt(rigamma(1,(K+nu)/2,(nu*lambda+sum(initak_prop^2))/2))
	initak_prop=rnorm(K)
	for (i in 1:K){
		sum_R_prop=sum(R_prop[groupings==gindex[i]])
		initak_prop[i]=initak_prop[i]*(bsig_prop*inittau_prop)/sqrt(nk*inittau_prop^2+bsig_prop^2)+(inittau_prop^2*sum_R_prop)/(nk*inittau_prop^2+bsig_prop^2)
	}
	tildey_prop=y-rep(initak_prop,each=nk)
	beta1=(2*base)/(max(tildey_prop)-min(tildey_prop))
	beta0=base-max(tildey_prop)*beta1
	transy=beta0+beta1*tildey_prop
	sigma=beta1*bsig_prop
	onedraw_prop=bart1draw(transy,onedraw_prop$mu_ij,sigma,sigma_mu,X,onedraw_prop$value,onedraw_prop$xindex,onedraw_prop$depth,onedraw_prop$varindex,onedraw_prop$variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
	gxprop=onedraw_prop$un_trans
	R_prop=y-gxprop
	psig_prop[post]=bsig_prop
	ptau_prop[post]=inittau_prop
	pgxprop[,post]=gxprop
	pakprop[,post]=initak_prop
}

#####Binary outcomes######
#Setup toy dataset
nk=5
K=50
num=nk*K
true_tau=1
p=10
X=matrix(runif(num*p),ncol = p)
Gx=1.35*(sin(pi*X[,1]*X[,2])+2*(X[,3]-.5)^2-X[,4]-0.5*X[,5])
A=rnorm(K,0,true_tau)
G_b=Gx+rep(A,each=nk)
Zik=rnorm(num)+G_b
pix=pnorm(G_b)
y=numeric(num)
for(i in 1:num) y[i]=ifelse(Zik[i]>0,1,0)
groupings=rep(1:K,each=nk)
gindex=sort(unique(groupings))

#Setup initial parameters
library(pscl)
library(zipfR)
library(truncnorm)
library(ROCR)
numpost=5000
numburnin=1000
initgx=as.numeric(apply(pgxreg,1,mean))
zlat=rtruncnorm(num,a=0,mean=initgx,sd=1)*as.numeric(y==1)+rtruncnorm(num,b=0,mean=initgx,sd=1)*as.numeric(y==0)
R_init=zlat-initgx
initak=numeric(K)
for(i in 1:K) initak[i]=mean(R_init[groupings==gindex[i]])
nu=2
lambda=1
alpha=0.95
beta=2
k=2
m=200
base=3
sigma_mu=base/(k*sqrt(m))
probgrow=0.5
probprune=0.5

#riBART
source("bart_fns.R")

#initializations
inittau_prop=sqrt(rigamma(1,(K+nu)/2,(nu*lambda+sum(initak^2))/2))
initak_prop=rnorm(K)

for (i in 1:K){
	sum_R_init=sum(R_init[groupings==gindex[i]])
	initak_prop[i]=initak_prop[i]*inittau_prop/sqrt(nk*inittau_prop^2+1)+(inittau_prop^2*sum_R_init)/(nk*inittau_prop^2+1)
}

tildey_prop=zlat-rep(initak_prop,each=nk)
beta1=(2*base)/(max(tildey_prop)-min(tildey_prop))
beta0=base-max(tildey_prop)*beta1
transy=beta0+beta1*tildey_prop
sigma=beta1
init_mu_ij=matrix(1/m*mean(transy),length(transy),m)
init_value=list()
init_xindex=list()
init_depth=list()
init_varindex=list()
init_variable=list()
for(i in 1:m){
	init_value[[i]]=NA
	init_xindex[[i]]=NA
	init_depth[[i]]=1
	init_varindex[[i]]=1
	init_variable[[i]]=rep(1,length(transy))	
}
onedraw_prop=bart1draw(transy,init_mu_ij,sigma,sigma_mu,X,init_value,init_xindex,init_depth,init_varindex,init_variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
gxprop=onedraw_prop$un_trans

#burnin
for(burn in 1:numburnin){
	tmpprop=gxprop+rep(initak_prop,each=nk)
	zlat_prop=rtruncnorm(num,a=0,mean=tmpprop,sd=1)*as.numeric(y==1)+rtruncnorm(num,b=0,mean=tmpprop,sd=1)*as.numeric(y==0)
	inittau_prop=sqrt(rigamma(1,(K+nu)/2,(nu*lambda+sum(initak_prop^2))/2))
	R_prop=zlat_prop-gxprop
	initak_prop=rnorm(K)
	for (i in 1:K){
		sum_R_prop=sum(R_prop[groupings==gindex[i]])
		initak_prop[i]=initak_prop[i]*inittau_prop/sqrt(nk*inittau_prop^2+1)+(inittau_prop^2*sum_R_prop)/(nk*inittau_prop^2+1)
	}
	tildey_prop=zlat_prop-rep(initak_prop,each=nk) 
	beta1=(2*base)/(max(tildey_prop)-min(tildey_prop))
	beta0=base-max(tildey_prop)*beta1
	transy=beta0+beta1*tildey_prop
	sigma=beta1
	onedraw_prop=bart1draw(transy,onedraw_prop$mu_ij,sigma,sigma_mu,X,onedraw_prop$value,onedraw_prop$xindex,onedraw_prop$depth,onedraw_prop$varindex,onedraw_prop$variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
	gxprop=onedraw_prop$un_trans
}

#posterior draws of tau, gx, and the random intercept
pgxprop=matrix(0,nrow=num,ncol=numpost)
pakprop=matrix(0,nrow=K,ncol=numpost)
ptau_prop=numeric(numpost)

for(post in 1:numpost){
	tmpprop=gxprop+rep(initak_prop,each=nk)
	zlat_prop=rtruncnorm(num,a=0,mean=tmpprop,sd=1)*as.numeric(y==1)+rtruncnorm(num,b=0,mean=tmpprop,sd=1)*as.numeric(y==0)
	inittau_prop=sqrt(rigamma(1,(K+nu)/2,(nu*lambda+sum(initak_prop^2))/2))
	R_prop=zlat_prop-gxprop
	initak_prop=rnorm(K)
	for (i in 1:K){
		sum_R_prop=sum(R_prop[groupings==gindex[i]])
		initak_prop[i]=initak_prop[i]*inittau_prop/sqrt(nk*inittau_prop^2+1)+(inittau_prop^2*sum_R_prop)/(nk*inittau_prop^2+1)
	}
	tildey_prop=zlat_prop-rep(initak_prop,each=nk) 
	beta1=(2*base)/(max(tildey_prop)-min(tildey_prop))
	beta0=base-max(tildey_prop)*beta1
	transy=beta0+beta1*tildey_prop
	sigma=beta1
	onedraw_prop=bart1draw(transy,onedraw_prop$mu_ij,sigma,sigma_mu,X,onedraw_prop$value,onedraw_prop$xindex,onedraw_prop$depth,onedraw_prop$varindex,onedraw_prop$variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
	gxprop=onedraw_prop$un_trans
	pgxprop[,post]=gxprop
	pakprop[,post]=initak_prop
	ptau_prop[post]=inittau_prop
}
