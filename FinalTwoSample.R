######### Required packages ########
require(clue, quietly=T)
require(energy, quietly = T)
require(randtoolbox, quietly = T)
require(pracma, quietly = T)
require(kernlab, quietly = T)
require(crossmatch, quietly = T)
require(HHG, quietly = T)
require(gTests, quietly = T)

######### Data generation ##########
m=200
n=200
data1=cbind(rcauchy(m,0,1),rcauchy(m,0,1))
data2=cbind(rcauchy(n,0.5,1),rcauchy(n,0,1))

######### Computing Rank Energy Statistic #########
computestatistic=function(x,y,m=nrow(x),n=nrow(y),dim=ncol(x),gridch=torus(m+n,dim))
{
  comdata=rbind(x,y)
  distmat=matrix(0,nrow=m+n,ncol=m+n)
  for(i in 1:(m+n))
    distmat[i,]=apply((comdata[i,]-t(gridch)),2,Norm)^2
  assignmentFUN=solve_LSAP(distmat)
  assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
  randenergySTAT=eqdist.etest(gridch[assignmentSOL[,2],],sizes = c(m,n), R=1)
  return(randenergySTAT$statistic)
}

######### Generating universal distribution ##########

gensamdist=function(M,N,dim,niter=30000,fixgrid=torus(M+N,dim))
{
  tstat=numeric(niter)
  for(i in 1:niter)
  {
    ranper=sample(M+N)
    print(i)
    tstat[i]=eqdist.etest(fixgrid[ranper,],sizes=c(M,N),R=1)$statistic
  }
  return(tstat)
}

########## Quantiles of null distribution ############
index=seq(0.5,1,0.01)
quanenergy1=quantile(mystat,index)

mystat=computestatistic(data1,data2)
tstat=gensamdist(m,n)
kmmd(data1,data2)$p.value

data=rbind(data1,data2)
rosdist=as.matrix(dist(data,diag=T,upper=T))
z=c(rep(0,m),rep(1,n))
crossmatchtest(z,rosdist)$pval
hhg.test.2.sample(rosdist,z,nr.perm=1000)$perm.pval.hhg.sl

mainf=function()
{
  perfvec=matrix(0,10,10)
  for(i in 1:10)
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
    #T2hhg=computestatisticrhhg(setdata1[,,i],setdata2[,,i])
    #perfvec[21,i]=length(which(quanenergyhhg>T2hhg))/100
    print(i)
  }
  for(i in 1:10)
  {
    for(j in 1:5)
    {
      perfvec[(3*j-1),i]=as.numeric(perfvec[(3*j),i]<0.05)
      perfvec[(3*j-2),i]=as.numeric(perfvec[(3*j),i]<0.1)
    }
    print(i)
  }
  return(perfvec)
}

#######################################


.libPaths("/rigel/stats/users/nd2560/rpackages")
library(data.table)
fileList=list.files("/rigel/home/nd2560",pattern="\\.RData$",full.name=T)
finout=array(0,dim=c(18,10,length(fileList)))
for(i in 1:length(fileList))
  finout[,,i]=get(load(fileList[i]))
meanmat=matrix(0,nrow=18,ncol=10)
for(i in 1:length(fileList))
  meanmat=meanmat+finout[,,i]
meanmat=meanmat/length(fileList)

sam=30
dat=matrix(0,sam,2)
mu1=c(-1,0)
mu2=c(3,0)
mu3=c(1,-1)
Sig1=matrix(c(1.25,-1,-1,1.25),2,2)/1.5
Sig2=matrix(c(1.25,1,1,1.25),2,2)/1.5
Sig3=matrix(c(1,0,0,0.25),2,2)/1.5

for(i in 1:sam)
{
  z=sample(3,1,prob=c(3/8,3/8,1/4))
  if(z==1)
    dat[i,]=mvrnorm(1,mu1,Sig1)
  if(z==2)
    dat[i,]=mvrnorm(1,mu2,Sig2)
  if(z==3)
    dat[i,]=mvrnorm(1,mu3,Sig3)
}
#sam=30
h=halton(sam,2)
dath=matrix(0,sam,2)
dath[,1]=dat[,1]/3
dath[,2]=dat[,2]/1.2
#dath[,2]=(dath[,2]+1.5)*2/3-0.5
#points(c(-0.2,-0.2,0.8,0.8),c(0.2,1.2,0.2,1.2))
#points(0,-1.4)
#points(h,pch=19,cex=0.8,col="blue")
#dat=dath
distmat=matrix(0,sam,sam)
for(i in 1:sam)
  distmat[i,]=apply((dath[i,]-t(h)),2,Norm)^2
assignmentFUN=solve_LSAP(distmat)
assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
nch=h[assignmentSOL[,2],]
nch[,1]=nch[,1]-0.2
nch[,2]=nch[,2]+0.3
plot(dath[,1]+0.1,dath[,2],pch=19,cex=1.15,col="red",xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(-1.25,1.45),xlim=c(-1,1.6))
rect(-0.2,0.3,0.8,1.3,border="green",lwd=3)
points(nch,pch=4,cex=1.4,col="blue")
arrows(dath[,1]+0.1,dath[,2],nch[,1],nch[,2],length=0,lty=4)
text(0.92,1.3,"(1,1)",cex=1.3)
text(-0.32,1.3,"(0,1)",cex=1.3)
text(-0.32,0.3,"(0,0)",cex=1.3)
text(0.92,0.3,"(1,0)",cex=1.3)
legend(-1.05,-0.4, c("Data points", "Empirical ranks"),bty='n',col = c("red", "blue", 6),pch = c(19, 4),cex=c(1.3,1.3))
datq=datq
nchq=nch

sam=30
x1=-runif(sam)
h=halton(sam,1)
plot(x1,pch=19,cex=0.8,col="red",xlab="X coordinate",ylab="Y coordinate",ylim=c(-1,1))
points(h,pch=19,cex=0.8,col="blue")
distmat=matrix(0,sam,sam)
for(i in 1:sam)
  distmat[i,]=(x1[i]-h)^2
assignmentFUN=solve_LSAP(distmat)
assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
nch=h[assignmentSOL[,2]]
Arrows(1:30,x1,1:30,nch,arr.adj=0.5)

n=10
h=(1:n)/n
x=rnorm(n,0,1)
plot(x,rep(2,n),pch=19,cex=1.2,col="red",xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(0.5,2.5),xlim=c(-1.5,1.7))
arrows(-1.5,2,1.7,2,length=0.1,code=3)
#abline(h=2,lty=1)
x=sort(x)
text(x[1],2.1,bquote("x"[(1)]),cex=1.3)
text(x[2],2.1,bquote("x"[(2)]),cex=1.3)
text(x[3],2.1,bquote("x"[(3)]),cex=1.3)
text(x[4],1.9,bquote("x"[(4)]),cex=1.3)
text(x[5],2.1,bquote("x"[(5)]),cex=1.3)
text(x[6],1.9,bquote("x"[(6)]),cex=1.3)
text(x[7],2.1,bquote("x"[(7)]),cex=1.3)
text(x[8]+0.07,1.9,bquote("x"[(8)]),cex=1.3)
text(x[9],2.1,bquote("x"[(9)]),cex=1.3)
text(x[10],2.1,bquote("x"[(10)]),cex=1.3)
#text(x[1],2.1,bquote("x"[(i)]),cex=0.6)
#text(sort(x),rep(2.1,n),labels=c(bquote("x"[(1)]),bquote("x"[(2)])),cex=0.8,font=0.5)
#points(h,rep(0,n),pch=19,cex=0.5,col="blue")
distmat=matrix(0,n,n)
for(i in 1:n)
  distmat[i,]=(x[i]-h)^2
assignmentFUN=solve_LSAP(distmat)
assignmentSOL=cbind(seq_along(assignmentFUN),assignmentFUN)
nch=h[assignmentSOL[,2]]
nch=nch*2.45-1.1
points(nch,rep(1,n),pch=4,cex=1.2,col="blue")
#abline(h=1,lty=2)
segments(-1.1,1,nch[10],1)
segments(-1.1,1.06,-1.1,0.94)
segments(-1.1,1.06,-1.07,1.06)
segments(-1.1,0.94,-1.07,0.94)
segments(nch[10],1.06,nch[10],0.94)
segments(nch[10],1.06,nch[10]-0.04,1.06)
segments(nch[10],0.94,nch[10]-0.04,0.94)
text(nch,rep(0.85,10),labels=h,cex=1.3)
text(-1.1,0.85,labels=0.0,cex=1.3)
arrows(x,rep(2,n),nch,rep(1,n),length=0,lty=3)
legend(-0.4,2.5, c("Data points"),bty='n',col = c("red"),pch = c(19),cex=c(1.3))
legend(-0.4,0.77, c("Empirical ranks"),bty='n',col = c("blue"),pch = c(4),cex=c(1.3))
#text(sort(x),rep(2.1,n),bquote("x"[(1)]),bquote("x"[(2)]),"x"[(3)],"x"[(4)],"x"[(5)],"x"[(6)],"x"[(7)],"x"[(8)],"x"[(9)],"x"[(10)]),cex=0.8,font=0.5)

d1=read.csv("checkquan-8-8.csv")[,-1]
d2=as.vector(read.csv("UniformInd-8-8.csv")[,-1])
dim(d1)
d=unlist(as.vector(d1[20,]))
quantile(d2,0.95)
qqplot(x=d,y=d2)
d1=as.matrix(d1,nrow=20,ncol=10000)
for(i in 1:29)
{
  v1[i]=ks.test(d1[i,],d1[30,])$p.val
}
v1=numeric(30)
for(i in 1:30)
{
  v1[i]=quantile(d1[i,],0.95)
}
s1=v1[c(2,6,10,14,18)]
s2=v1[c(2,6,10,14,18)]
s3=v1[c(2,6,10,14,18)]
s4=v1[c(2,6,10,14,18)]
s5=v1[c(2,6,10,14,18)]
s6=v1[c(2,6,10,14,18)]
qqplot(d1[12,],d1[20,],pch=16,cex=0.5)
lines(seq(0.1,1.6,0.01),seq(0.1,1.6,0.01),ty='l',lwd=2,col="red")
