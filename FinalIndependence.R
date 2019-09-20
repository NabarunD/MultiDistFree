######### Required packages ########
require(clue, quietly=T)
require(energy, quietly = T)
require(randtoolbox, quietly = T)
require(pracma, quietly = T)
require(HHG, quietly = T)
require(dHSIC, quietly = T)
#require(copula, quietly = T)
require(IndepTest, quietly = T)
require(SpatialNP, quietly = T)
require(EnvStats, quietly = T)

######### Data generation ##########

quanenergydcov=read.csv("QuanenergyRankDcov-100-5-5.csv")[,-1]
quanenergyhhg=read.csv("QuanenergyRankHHG-100-5-5.csv")[,-1]
n=100
d=1
s=18
setdata1=array(0,dim=c(n,d,s))
setdata2=array(0,dim=c(n,d,s))

#### Cauchy setup

for(i in 1:d)
{
prep=rnorm(n,0,1)
setdata1[,i,1]=0.2*rcauchy(n,0,1)+prep
setdata2[,i,1]=prep+0.2*rcauchy(n,0,1)
}

#### Parabola

for(i in 1:d)
{
setdata1[,i,2]=runif(n,-1,1)
setdata2[,i,2]=(setdata1[,i,2]^2+runif(n,0,1))/2
}

#### Circle

for(i in 1:d)
{
prep=runif(n,-1,1)
setdata1[,i,3]=sin(pi*prep)+rnorm(n,0,1)/8
setdata2[,i,3]=cos(pi*prep)+rnorm(n,0,1)/8
}

#### Two parabolas

for(i in 1:d)
{
setdata1[,i,4]=runif(n,-1,1)
A=rbinom(n,1,prob=0.5)
setdata2[,i,4]=A*(setdata1[,i,4]^2+runif(n,0,1))/2+(1-A)*(setdata1[,i,4]^2+runif(n,0,1))*(-0.5)
}

#### Diamond

theta=pi/4
for(i in 1:d)
{
prep1=runif(n,-1,1)
prep2=runif(n,-1,1)
d1=sin(theta)*prep1+cos(theta)*prep2
d2=-sin(theta)*prep1+cos(theta)*prep2
#setdata1[,i,5]=sin(theta)*prep1+cos(theta)*prep2
#setdata2[,i,5]=-sin(theta)*prep1+cos(theta)*prep2
}

#### W-shape

d1=numeric(n)
d2=numeric(n)
for(i in 1:d)
{
prep=runif(n,-1,1)
prep1=runif(n,0,1)
prep2=runif(n,0,1)
d1=prep+prep1/3
d2=4*(prep^2-0.5)^2+prep2
#setdata1[,i,6]=prep+prep1/3
#setdata2[,i,6]=4*(prep^2-0.5)^2+prep2
}

#### Gaussian mixture

for(i in 1:d)
{
dx <- rnorm(n)/3
dy <- rnorm(n)/3
cx <- sample( c(-1,1), size=n, replace=T )
cy <- sample( c(-1,1), size=n, replace=T )
setdata1[,i,7] <- cx + dx
setdata2[,i,7] <- cy + dy 
}

#### Mutual Information 1

l = 1
for(j in 1:d)
{
i = 1
while(i<=n)
{
  prep=runif(2,-pi,pi)
  A=rbinom(1,1,prob=0.5*(1+sin(l*prep[1])*sin(l*prep[2])))
  if(A==1)
  {
    setdata1[i,j,8]=prep[1]
    setdata2[i,j,8]=prep[2]
    i=i+1
  }
}
}

#### Mutual Information 2

l = 5
for(i in 1:d)
{
L = sample(l,n,replace=T)
Theta = runif(n,0,2*pi)
ep1 = rnorm(n,0,1)
ep2 = rnorm(n,0,1)
setdata1[,i,9]=L*cos(Theta)+ep1/4
setdata2[,i,9]=L*sin(Theta)+ep2/4
}

#### Mutual Information 3

for(i in 1:d)
{
setdata1[,i,10]=runif(n,-1,1)
ep=rnorm(n,0,1)
rh=2
setdata2[,i,10]=ep*(abs(setdata1[,i,10]))^(rh)
}

##### Linear relationship

for(i in 1:d)
{
  prep=rnorm(n,0,1)
  setdata1[,i,11]=runif(n,0,1)
  setdata2[,i,11]=setdata1[,i,11]+prep
}

##### Fractional Exponent

for(i in 1:d)
{
  prep=rnorm(n,0,1)
  setdata1[,i,12]=runif(n,0,1)
  setdata2[,i,12]=setdata1[,i,12]^(1/4)+prep
}

#### Cauchy setup independent

for(i in 1:d)
{
  #prep=rnorm(n,0,1)
  setdata1[,i,13]=rcauchy(n,0,0.5)
  setdata2[,i,13]=rcauchy(n,0,0.5)
}

#### Pareto setup independent

for(i in 1:d)
{
  prep=rnorm(n,0,1)
  setdata1[,i,14]=prep+rpareto(n,1,2)^2
  setdata2[,i,14]=prep+rpareto(n,1,2)^2
}

#### Logarithmic setup

for(i in 1:d)
{
  #prep=rnorm(n,0,1)
  setdata1[,i,15]=rnorm(n,0,1)
  setdata2[,i,15]=log(setdata1[,i,15]^2)
}

#### Epsilon setup

for(i in 1:d)
{
  #prep=rnorm(n,0,1)
  prep=rnorm(n,0,1)
  setdata1[,i,16]=rnorm(n,0,1)
  setdata2[,i,16]=prep*setdata1[,i,16]
}

#### Quadratic setup

for(i in 1:d)
{
  #prep=rnorm(n,0,1)
  prep=rnorm(n,0,3)
  setdata1[,i,17]=rnorm(n,0,1)
  if(i<=2)
  setdata2[,i,17]=setdata1[,i,17]+4*setdata1[,i,17]^2+prep
  if(i>2)
  setdata2[,i,17]=prep
}

#### Multivariate Normal setup

mn=rep(0,2*d)
Sig=diag(2*d)
for(i in 1:d)
{
  for(j in (d+1):(2*d))
  {
    Sig[i,j]=0.1
    Sig[j,i]=0.1
  }
}
totdata=mvrnorm(n,mn,Sig)
setdata1[,,18]=totdata[,(1:d)]
setdata2[,,18]=totdata[,((d+1):(2*d))]

######### Computing Rank Energy Statistic #########
computestatisticrdcov=function(x,y,dim1=ncol(x),dim2=ncol(y),n=nrow(x),gridch=halton(n,dim1+dim2),gridch1=as.matrix(gridch[,(1:dim1)]),gridch2=as.matrix(gridch[,((dim1+1):(dim1+dim2))]))
{
  distmat1=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat1[i,]=apply((x[i,]-t(gridch1)),2,Norm)^2
  assignmentFUN1=solve_LSAP(distmat1)
  assignmentSOL1=cbind(seq_along(assignmentFUN1),assignmentFUN1)
  distmat2=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat2[i,]=apply((y[i,]-t(gridch2)),2,Norm)^2
  assignmentFUN2=solve_LSAP(distmat2)
  assignmentSOL2=cbind(seq_along(assignmentFUN2),assignmentFUN2)
  randcovSTAT=dcov.test(gridch1[assignmentSOL1[,2],],gridch2[assignmentSOL2[,2],], R=1)
  return(randcovSTAT$statistic)
}

######### Generating universal distribution ##########
gensamdistrdcov=function(N,dim1,dim2,niter=15000,fixgrid=halton(N,dim1+dim2),fixgrid1=as.matrix(fixgrid[,(1:dim1)]),fixgrid2=as.matrix(fixgrid[,((dim1+1):(dim1+dim2))]))
{
  tstat=numeric(niter)
  for(i in 1:niter)
  {
    ranper=sample(N)
    tstat[i]=dcov.test(fixgrid1,fixgrid2[ranper,],R=1)$statistic
    print(i)
  }
  return(tstat)
}

######### Computing Rank HHG Statistic #########
computestatisticrhhg=function(x,y,dim1=ncol(x),dim2=ncol(y),gridch=halton(n,dim1+dim2),gridch1=as.matrix(gridch[,(1:dim1)]),gridch2=as.matrix(gridch[,((dim1+1):(dim1+dim2))]))
{
  distmat1=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat1[i,]=apply((x[i,]-t(gridch1)),2,Norm)^2
  assignmentFUN1=solve_LSAP(distmat1)
  assignmentSOL1=cbind(seq_along(assignmentFUN1),assignmentFUN1)
  distmat2=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat2[i,]=apply((y[i,]-t(gridch2)),2,Norm)^2
  assignmentFUN2=solve_LSAP(distmat2)
  assignmentSOL2=cbind(seq_along(assignmentFUN2),assignmentFUN2)
  Dx = as.matrix(dist(gridch1[assignmentSOL1[,2],]), diag = TRUE, upper = TRUE)
  Dy = as.matrix(dist(gridch2[assignmentSOL2[,2],]), diag = TRUE, upper = TRUE)
  randcovSTAT=hhg.test(Dx,Dy, nr.perm=1)$sum.lr/(n^2)
  return(randcovSTAT)
}


######### Generating universal distribution ##########
gensamdistrhhg=function(N,dim1,dim2,niter=5000,fixgrid=halton(N,dim1+dim2),fixgrid1=as.matrix(fixgrid[,(1:dim1)]),fixgrid2=as.matrix(fixgrid[,((dim1+1):(dim1+dim2))]))
{
  tstat=numeric(niter)
  for(i in 1:niter)
  {
    ranper=sample(N)
    Dx = as.matrix(dist(fixgrid1), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist(fixgrid2[ranper,]), diag = TRUE, upper = TRUE)
    tstat[i]=hhg.test(Dx,Dy,nr.perm=1)$sum.lr/(N^2)
    print(i)
  }
  return(tstat)
}

######### Quantiles of null distribution ############
#index=seq(0.01,1,0.01)
#quanenergy2=quantile(mystat8,index)

######### Other tests ##########
#sr.indep.test(data1,data2,score = c("symmsign"))
#sr.indep.test(data1,data2,score = c("sign"))
#sr.indep.test(data1,data2,score = c("rank"))
#dcov.test(data1,data2,R=200)
#dhsic.test(data1,data2)
#MINTauto(data1,data2,kmax=10)
#MINTav(data1,data2,K=1:99)
#MINTperm(data1,data2,k=99)
#Dx = as.matrix(dist(set9data1), diag = TRUE, upper = TRUE)
#Dy = as.matrix(dist(set9data2), diag = TRUE, upper = TRUE)
#hhg.test(Dx,Dy,nr.perm = 1000)$perm.pval.hhg.sl
#length(which(mystat2>computestatistic(as.matrix(data1),as.matrix(data2))))/30000
mainf=function()
{
  perfvec=matrix(0,21,18)
  for(i in 1:18)
  {
    perfvec[3,i]=as.numeric(sr.indep.test(setdata1[,,i],setdata2[,,i],score=c("rank"))$p.val)
    perfvec[6,i]=as.numeric(dcov.test(setdata1[,,i],setdata2[,,i],R=500)$p.val)
    perfvec[9,i]=as.numeric(dhsic.test(setdata1[,,i],setdata2[,,i])$p.val)
    perfvec[12,i]=MINTav(setdata1[,,i],setdata2[,,i],K=1:(n-1))
    Dx = as.matrix(dist(setdata1[,,i]), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist(setdata2[,,i]), diag = TRUE, upper = TRUE)
    perfvec[15,i]=as.numeric(hhg.test(Dx,Dy,nr.perm = 1000)$perm.pval.hhg.sl)
    T1dcov=computestatisticrdcov(setdata1[,,i],setdata2[,,i])
    perfvec[18,i]=length(which(quanenergydcov>T1dcov))/100
    T2hhg=computestatisticrhhg(setdata1[,,i],setdata2[,,i])
    perfvec[21,i]=length(which(quanenergyhhg>T2hhg))/100
    print(i)
  }
  for(i in 1:18)
  {
    for(j in 1:7)
    {
      perfvec[(3*j-1),i]=as.numeric(perfvec[(3*j),i]<0.05)
      perfvec[(3*j-2),i]=as.numeric(perfvec[(3*j),i]<0.1)
    }
    print(i)
  }
libPaths("/rigel/stats/users/nd2560/rpackages")
library(data.table)
fileList=list.files("/rigel/home/nd2560",pattern="\\.RData$",full.name=T)
finout=array(0,dim=c(30,10,length(fileList)))
for(i in 1:length(fileList))
  finout[,,i]=get(load(fileList[i]))
meanmat=matrix(0,nrow=20,ncol=10)
for(i in 1:length(fileList))
  meanmat=meanmat+finout[,,i]
meanmat=meanmat/length(fileList)
}

######## Combine Output ############



.libPaths("/rigel/stats/users/nd2560/rpackages")
library(data.table)
fileList=list.files("/rigel/home/nd2560",pattern="\\.RData$",full.name=T)
finout=NULL
for(i in 1:length(fileList))
  finout=c(finout,get(load(fileList[i])))
#setwd("/rigel/stats/users/nd2560")

write.csv(finout,"OutputFile.csv")


##########################################
n=200
d=3
#iter=100
#set2check=numeric(iter)
#for(i in 1:iter)
#{

for(i in 1:d)
{
  prep=runif(n,-1,1)
  setdata1[,i,3]=sin(pi*prep)+rnorm(n,0,1)/8
  setdata2[,i,3]=cos(pi*prep)+rnorm(n,0,1)/8
}

#  T1dcov=computestatisticrdcov(setdata1[,,2],setdata2[,,2])
#  set2check[i]=length(which(quanenergydcov>T1dcov))/100
#  print(i)
#  print(T1dcov)
#}

t1=rcauchy(200,0,1)
t2=rcauchy(200,0,1)
A=rbinom(200,1,0.05)
sum(A)
t2=A*t1+(1-A)*t2
#dcov.test(t1,t2,R=200)
dhsic.test(t1,t2)$p.value


####################################################
dim=3
d1=matrix(0,200,3)
d2=matrix(0,200,3)
for(i in 1:dim)
{
  d1[,i]=rnorm(200,0,1)
  d2[,i]=rnorm(200,0,1)
}
A=rbinom(200,1,0.04)
d2=(1-A)*d2+A*d1
dcov.test(d1,d2,R=200)$p.val
dhsic.test(d1,d2)$p.val
length(which(quanenergy2>computestatisticrdcov(d1,d2)))/100

######################################################

tab=read.csv("Independence-200-3-3-revision-5.csv")[,-1]
tab=t(tab)
tab1=tab[,-c(3,6,9,12,15,18)]
print(xtable(tab1))

n=200
x1=rnorm(n,0,1)
d1=abs(x1+rpareto(n,1,1))^{1.5}
d2=abs(x1+rpareto(n,1,1))^{1.5}

n=200
x1=rnorm(n,0,1)
d1=abs(x1+rpareto(n,1,2)^2)^{1}
d2=abs(x1+rpareto(n,1,2)^2)^{1}

d1=as.vector(setdata1[,,6])
d2=as.vector(setdata2[,,6])
hseq=halton(n,2)
hseq1=hseq[,1]
hseq2=hseq[,2]
distmat1=matrix(0,nrow=n,ncol=n)
for(i in 1:n)
  distmat1[i,]=sapply((d1[i]-hseq1),Norm)^2
assignmentFUN1=solve_LSAP(distmat1)
assignmentSOL1=cbind(seq_along(assignmentFUN1),assignmentFUN1)
distmat2=matrix(0,nrow=n,ncol=n)
for(i in 1:n)
  distmat2[i,]=sapply(d2[i]-hseq2,Norm)^2
assignmentFUN2=solve_LSAP(distmat2)
assignmentSOL2=cbind(seq_along(assignmentFUN2),assignmentFUN2)
plot(hseq1[assignmentSOL1[,2]],hseq2[assignmentSOL2[,2]])
plot(d1,d2,xlim=c(0,15),ylim=c(0,15))
kk=ksmoo

x=read.csv("GaussianIndependence-200-3-3-3.csv")[,-1]
#y1=x[c(20:11,1:10),c(1,2,9,10,3,4,5,6)]
y1=x[,c(4,5,7,8,10,11,13,14,16,17)]
y1=rbind(y1[20:11,],0,y1[1:10,])
rh=c(0.075,0.15,0.225,0.3,0.375,0.450,0.525,0.6,0.75,0.9)
rh=c(rh,-rh)
rh=c(rh[20:11],0,rh[1:10])
newdata=seq(-0.9,0.9,0.025)
l1=loess(y1[,2]~rh,data=y1,span=0.5)
sml1=predict(l1,newdata)
l2=loess(y1[,2]~rh,data=y1,span=0.5)
sml2=predict(l2,newdata)
l3=loess(y1[,4]~rh,data=y1,span=0.5)
sml3=predict(l3,newdata)
l4=loess(y1[,6]~rh,data=y1,span=0.5)
sml4=predict(l4,newdata)
l5=loess(y1[,8]~rh,data=y1,span=0.5)
sml5=predict(l5,newdata)
l6=loess(y1[,10]~rh,data=y1,span=0.5)
sml6=predict(l6,newdata)
sml2[which(sml2>1)]=1
sml3[which(sml3>1)]=1
sml4[which(sml4>1)]=1
sml5[which(sml5>1)]=1
sml6[which(sml6>1)]=1
par(bg="white")
plot(newdata,sml2,ty='l',lwd=2,xlab="Correlation",ylab="Power",cex.lab=1.5,xlim=c(-0.75,0.75),ylim=c(0,1.1),col="yellow")
lines(newdata,sml2,ty='l',lwd=2,xlab="Correlation",ylab="Power",cex.lab=1.5,col="red")
lines(newdata,sml3,ty='l',lwd=1,xlab="Correlation",ylab="Power",cex.lab=1.5,col="blue")
lines(newdata,sml4,ty='l',lwd=2,xlab="Correlation",ylab="Power",cex.lab=1.5,col="green")
lines(newdata,sml5,ty='l',lwd=2,xlab="Correlation",ylab="Power",cex.lab=1.5,col="violet")
lines(newdata,sml6,ty='l',lwd=2,xlab="Correlation",ylab="Power",cex.lab=1.5,col="black")
legend("bottomright",c("DCoV","HSIC","MINT","HHG","RdCov"),col=c("red","blue","green","violet","black"),lwd=2,cex=1.5,lty=1)
View(y1)

x=read.csv("uniformpred-2-2-1.csv")[,-1]
y=read.csv("checkquan-2-2.csv")[,-1]
dim(y)
z=unlist(y[20,])
plot(density(x),lwd=2,col="blue",xlab="",cex.lab=1.5,main="")
lines(density(z),lwd=2,col="red",xlab="",cex.lab=1.5,main="")
legend("topright",c("n*DCoV^2 for n=1500 (uniform)","n*RdCov^2 for n=1000"),col=c("blue","red"),cex=1.1,lwd=2,lty=1)
x=read.csv("UniformInd-8-8.csv")[,-1]
y=read.csv("checkquan-8-8.csv")[,-1]
x=read.csv("Twosam-1-8-quantile.csv")[,-1]
x=read.csv("UniformTwoSam-2-2-1.csv")[,-1]
y=read.csv("twosam-2-2-1.csv")[,-1]
plot(density(x),lwd=2,col="blue",xlab="",cex.lab=1.5,main="")
lines(density(z),lwd=2,col="red",xlab="",cex.lab=1.5,main="")
legend("topright",c("Scaled EN^2 for m=n=1500 (uniform)","Scaled RE^2 for n=1000"),col=c("blue","red"),bty="n",cex=1.2,lwd=2,lty=1)
x=read.csv("twosam-8-8-1.csv")[,-1]
z=unlist(x[20,])
x=read.csv("UniformTwoSam-8.csv")[,-1]
y=read.csv("twosam-8-8-1.csv")[,-1]
