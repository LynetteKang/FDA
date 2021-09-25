###############################  Project 1 ####################################

## load dataset: fly.Rdata
load("C:/Users/lenovo/Desktop/2021_函数型数据/Data sets/fly.RData")
y.egg = medfly$eggcount # eggcount matrix among 26 days for 50 flies
y.life = medfly$lifetime # lifetime of 50 flies


## About dataset:
#   The Fruit lies data consist of records of the number of eggs laid by 50 fruit
#   fies on several days, along with each individual's total lifespan.


##  50只果蝇生存时间直方图
#

par(mfrow=c(1,1))
plot(x = 1:length(y.life), y = y.life, xlab = "Index of 50 Flies", 
     ylab = "Life span", cex.lab=1.2,cex.axis=1.5,main = "Life time of fruit flies")
#lines(x = 1:length(y.life), y =y.life.hat,lwd=2,col="black") # fitting line
abline(h=c(mean(y.life),min(y.life),max(y.life)),lty=2,lwd=2,col=c("red","gray","gray"))# average line

which.max(y.life) # N0.38 has max life time = 739, died late
which.min(y.life) # N0.49 has min life time = 124, died early

## 寿命长短分布情况
# ggplot
ggplot(lifetime, aes(x=longitivity)) + 
  geom_histogram(aes(y=..density..),     
                 binwidth=20,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+ # 概率分布y轴density
  ggtitle("Histogram of Lifespan after laying eggs")+
  xlab("Longitivity/day")+ #加标题坐标
  ylab("Density")+
  theme( # 设置标题格式
    plot.title = element_text(color="Black", size=16, face="bold.italic"),
    axis.title.x = element_text(color="#993333", size=14, face="bold"),
    axis.title.y = element_text(color="#993333", size=14, face="bold")
  ) 
  

## 画25条时间和产蛋量的图
#
library("ggplot2")
library("reshape2")
egg <- data.frame(y.egg)
egg <-cbind(1:26,egg)
colnames(egg)<-c("Day", sapply(seq(1, 50), function(x){paste("Fly", x, sep="")}))
egg <- melt(egg,id.vars="Day",value.name = "value",variable.name = "variable")
p <-ggplot(data=egg,aes(x=Day,y=value,group=variable))+
  geom_line(size=0.1,alpha=0.5)+
  ggtitle("Eggcount with days for 50 Fruit Flies")+
  xlab("Day") + ylab("Egg count")+
  theme( # 设置标题格式
  plot.title = element_text(color="Black", size=16, face="bold.italic"),
  axis.title.x = element_text(color="#993333", size=14, face="bold"),
  axis.title.y = element_text(color="#993333", size=14, face="bold")
) 
#添加均值线
p<- p+stat_summary(aes(x=Day, y=value, group=1), 
               fun=mean, colour='red', geom = 'line',size=1.2)
#添加第一个果蝇的变化图
p<- p+geom_line(data=subset(egg, variable=="Fly1"), 
                aes(x=Day, y=value), colour="blue", size=1.2)
p
#添加图例
p + scale_colour_manual("",
                        breaks = c("Mean Fly community","50 Fruit Flies","First Fruit Fly"),
                        values = c("red","grey","blue")) +
  theme(legend.title=element_blank(),
        legend.position="top")
  
#x y axis坐标范围




## 50只果蝇分布egg count & day
# 此处y.egg 是matrix
par(mfrow=c(3,5))
for (i in 46:50) {
  plot(x = 1:26, y =y.egg[,i], xlab = "Day", 
       ylab = "Egg count", cex.lab=1.2,cex.axis=1.5,main = paste(i,"th Fly"))
  lines(x = 1:26, y =y.egg[,i])
}


## Non-parametric Smoothing Method:
## Set up true value and X(t) basis matrix
y.egg = medfly$eggcount
y.fly1 = y.egg[,1] # number 1 fly
y= apply(y.egg, 1, mean) # use median to represent egg count
X = matrix(1:26,26,10)

#### 1.Plynomial basis ####
## Define basis funtion:

for(i in 1:10){ X[,i] = ((X[,i]-14)/8)^(i-1)}  # kind of standarlization
head(X)
X[,1] = 0.2

par(mfrow=c(2,5))
for(i in 1:10){
  yhat = X[,1:i]%*%solve(t(X[,1:i])%*%X[,1:i])%*%(t(X[,1:i])%*%y)
  #matplot(150:216,X[150:216,1:i],type='l',lwd=2,lty=1,xlab='day',ylab='basis',cex.lab=1.5,cex.axis=1.5)
  plot(1:26,y,col=2,cex=1.5,xlab='day',ylab='Eggcount',cex.lab=1.5,
       main=paste(i,"Polynomial basis"),cex.axis=1.5)
  lines(1:26,yhat,lwd=2,col=4)
}


#### 2.Fourier basis ####

## Define basis funtion:
daybasis26 <- create.fourier.basis(rangeval=c(0, 26), nbasis = 26)
bvals = eval.basis(1:26,daybasis26) 

## evaluate \Phi_i(t) for each t
par(mfrow=c(1,5),ask=F)
for(i in 1:5){
  X = bvals[,1:(2*i+1)]  
  yhat = X%*%solve(t(X)%*%X)%*%(t(X)%*%y)
  #matplot(1:26,X[,(2*i):(2*i+1)],type='l',lwd=2,lty=1,xlab='day',ylab='basis',cex.lab=1.5,cex.axis=1.5)
  plot(1:26,y,col=2,cex=1.5,xlab='day',ylab='Eggcount',cex.lab=1.5,
       main=paste(2*i+1,"Fourier Basis"),cex.axis=1.5)
  lines(1:26,yhat,lwd=2,col=4)
}

#### 3.B-spline ####

## knots
# spline basis with 11 interior knots
# 12 months, 13 knots 
# total number of the basis functions is norder+ #interior knots = n + 11.
knots = c(seq(0,26,2))
bbasis = list() 

help(create.bspline.basis)
bbasis[[1]] = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=1)
bbasis[[2]] = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=2)
bbasis[[3]] = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=3)
bbasis[[4]] = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=4) #The default of 4 gives cubic splines.
bbasis[[5]] = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=5)
bbasis[[6]] = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=6)

# plot the data with fitted curves
# fda function: data2fd(y,1:365,bbasis[[5]])
#               to calculate the fitted curve for y,  
# or use the following direct way:  
#               (1)calculate the estimate the coefficients of the spline basis functions 
#               (2)calculate the fitted values.

par(mfrow=c(2,3))
for(i in 1:6){
  plot(1:26,y,col=2,ylab='Eggcount',xlab='day',
       main=paste(i,"order B-Spline"),cex.lab=1.5,cex.axis=1.5) # true values
  abline(v=knots,lty=2,lwd=1)
  lines(Data2fd(y,argvals=1:26,basisobj=bbasis[[i]]),col=4,lwd=2) # fitted values
  #bvals = eval.basis(1:365,bbasis[[i]]) 
  #plot(1:365,bvals[,i+1],type='l',lty=1,lwd=2,col=2,xlab='day',ylab='basis',cex.lab=1.5,cex.axis=1.5)
  #abline(v=knots,lty=2,lwd=1)
}


#### 4.Penalty: D2-operator ####

library(fda)
knots = c(seq(0,26,3),26) # knots each 3 days
bbasis = create.bspline.basis(rangeval=c(0,26),breaks=knots,norder=6) # b-spline with order 6
y.egg = medfly$eggcount
y.prec= apply(y.egg, 1, mean) # replace VancPrec
day = 1:26


# least-squares smooth:
precfd = smooth.basis(day,y.prec,bbasis)
plotfit.fd(y.prec,day,precfd$fd)

## define a linear differential operator (Lfd) object
# simplest Lfds: penalize the second derivative x (acceleration) 二阶导惩罚
D2lfd = int2Lfd(2) # 定义二阶求导算子
D2lfd
names(D2lfd)

#evaluate an Lfd
plot(precfd$fd,Lfdobj=D2lfd) # 求导之后的图像


# plotting the functional derivative
plot(deriv.fd(precfd$fd,2))


## define general notions of smoothing 
# vec2Lfd(): Make a Linear Differential Operator Object from a Vector
harmLfd = vec2Lfd(c(0,(2*pi/26)^2,0), c(0, 26)) # L operator 
harmLfd

# 思考题1 不管了


# 
plot(precfd$fd,Lfdobj=harmLfd)

#
eggfd = smooth.basis(day,y.egg,bbasis)
#
plot(eggfd)
plot(eggfd,Lfdobj=harmLfd)

# ends tend to be rough, examine the middle of a function
plot(eggfd,Lfdobj=harmLfd,y=-5000:5000)
abline(v=knots,col=2,lty=2)


## Saturated Bspline basis 饱和模型后续做gcv选lambda 和K

# bspline 建立saturated model with order 6
fbasis = create.bspline.basis(c(0,26),nbasis=30,norder=6) # nbasis = ninterior + order


## Define an fdPar (functionalparameter) object.
par(mfrow = c(2,1))
D2fdPar = fdPar(fbasis,Lfdobj=int2Lfd(2),lambda=1e4)  # lambda infinity
#验证了当lambda 无穷大的时候，线性模型
sprecfd = smooth.basis(day,y.prec,D2fdPar) #replace the basis object in smooth.basis with D2fdPar:
plotfit.fd(y.prec,day,sprecfd$fd,main= expression("penalty:" ~ L ==  omega^2 * D + D^3),ylab = "Eggcount", xlab="Day") 
text(0, 25, adj=0, 'lambda tends to infinity',col = "blue")


# 思考题2 重新做如下定义，Does the above work? How can you rectify this if not?
harmfdPar= fdPar(fbasis,Lfdobj=harmLfd,lambda=1e4)  # 定义L operation (harmonic)
sharmprecfd = smooth.basis(day,y.prec,harmfdPar)
#Compare the results by looking at
plotfit.fd(y.prec,day,sharmprecfd$fd,main= expression("penalty:" ~ L ==  omega^2 * D + D^3),ylab = "Eggcount", xlab="Day")
text(0, 25, adj=0, 'lambda tends to infinity',col = "blue")
#plot(sharmprecfd$fd,add=T,col=2)



## D2: change with lambda

gcv = rep(0,8)
df = rep(0,8)
sse = rep(0,8)

par(mfrow = c(2,4))
for(i in 1:8){
  lambda=10^{i-4}
  D2fdPar = fdPar(fbasis,Lfdobj=int2Lfd(2),lambda=lambda)
  sprecfd = smooth.basis(day,y.prec,D2fdPar)
  gcv[i] = sprecfd$gcv # gcv estimation
  df[i] = sprecfd$df   # df estimation 
  sse[i] = sprecfd$SSE  # sse estimation
  mainstr = paste('lambda = ',lambda,' df = ',round(df[i],2),' gcv = ',round(gcv[i]),sep='') 
  par(ask=TRUE)
  plot(day,y.prec,col=2,xlab='day',ylab='Egg Count',
       main=mainstr,cex.lab=1.5,cex.axis=1.5)
  lines(sprecfd$fd,col=4,lwd=2)
}

par(mfrow = c(1,3))
plot(-3:4,gcv,type="b",xlab = expression(log(lambda)),col="red",cex.lab=1.5,cex.axis=1.5,lwd=2)
abline(v=log10(10^(which.min(gcv)-4)),lty=2,lwd=2) # log10(lambda) =1,lambda = 10 
text(1.3,160, adj=0, 'min(GCV): 11.26',cex=1.2)
text(1.3,150, adj=0, expression(lambda == 10),cex=1.2)

plot(-3:4,sse,type="b", xlab = expression(log(lambda)),col="dark green",cex.lab=1.5,cex.axis=1.5,lwd=2)
abline(h = sse[which.min(gcv)],v=1,lty=2,lwd=2) # log10(lambda) =1,sse=170.9366
abline(v=1,lty=2,lwd=2)
text(-2.5,750, adj=0, 'sse: 170.9366',cex=1.2)
text(-2.5,700, adj=0, expression(lambda == 10),cex=1.2)

plot(-3:4,df,type="b", xlab = expression(log(lambda)),col="blue",cex.lab=1.5,cex.axis=1.5,lwd=2)
abline(h = df[which.min(gcv)],v=1,lty=2,lwd=2) # log10(lambda) =1,df=6.131637
abline(v=1,lty=2,lwd=2) 
text(2,22.5, adj=0, 'df: 6.1316',cex=1.2)
text(2,21, adj=0, expression(lambda == 10),cex=1.2)


#### 5.Penalty: L-operator ####

## L: change with lambda

# Saturated model for L:
harmfdPar= fdPar(fbasis,Lfdobj=harmLfd,lambda=1e4)  # 定义L operation (harmonic)
sharmprecfd = smooth.basis(day,y.prec,harmfdPar)

gcv = rep(0,8)
df = rep(0,8)
sse = rep(0,8)

par(mfrow = c(2,4))
for(i in 1:8){
  lambda=10^{i-4}
  harmfdPar = fdPar(fbasis,Lfdobj=harmLfd,lambda=lambda)
  sharmprecfd = smooth.basis(day,y.prec,harmfdPar)
  gcv[i] = sharmprecfd$gcv # gcv estimation
  df[i] = sharmprecfd$df   # df estimation 
  sse[i] = sharmprecfd$SSE  # sse estimation
  mainstr = paste('lambda = ',lambda,' df = ',round(df[i],2),' gcv = ',round(gcv[i]),sep='') 
  par(ask=TRUE)
  plot(day,y.prec,col=2,xlab='day',ylab='Egg Count',
       main=mainstr,cex.lab=1.5,cex.axis=1.5)
  lines(sharmprecfd$fd,col=4,lwd=2)
}

par(mfrow = c(1,3))
plot(-3:4,gcv,type="b",xlab = expression(log(lambda)),col="red",cex.lab=1.5,cex.axis=1.5,lwd=2)
abline(v=log10(10^(which.min(gcv)-4)),lty=2,lwd=2) # log10(lambda) =2,lambda = 100 
text(-2,60, adj=0, 'min(GCV): 10.32',cex=1.2)
text(-2,57, adj=0, expression(lambda == 100),cex=1.2)

plot(-3:4,sse,type="b", xlab = expression(log(lambda)),col="dark green",cex.lab=1.5,cex.axis=1.5,lwd=2)
abline(h = sse[which.min(gcv)],v=1,lty=2,lwd=2) # log10(lambda) =1,sse=163.0236
abline(v=log10(10^(which.min(gcv)-4)),lty=2,lwd=2)
text(-2,430, adj=0, 'sse: 163.02',cex=1.2)
text(-2,400, adj=0, expression(lambda == 100),cex=1.2)

plot(-3:4,df,type="b", xlab = expression(log(lambda)),col="blue",cex.lab=1.5,cex.axis=1.5,lwd=2)
abline(h = df[which.min(gcv)],v=1,lty=2,lwd=2) # log10(lambda) =1,df=5.738467
abline(v=log10(10^(which.min(gcv)-4)),lty=2,lwd=2) 
text(2.4,22, adj=0, 'df: 5.74',cex=1.2)
text(2.4,21, adj=0, expression(lambda == 100),cex=1.2)


## Part 3: simulation

#### 1. Choose the number of basis functions ####
# focus on the daily precipitation at Vancouver. 
# use Fourier basis to fit the data and its second derivative using the different number of basis functions.

library('fda')
y = CanadianWeather$dailyAv[,'Vancouver','Precipitation.mm']
daybasis365 <- create.fourier.basis(c(0, 365), 365)
par(mfrow=c(1,2),ask=T)
bvals = eval.basis(1:365,daybasis365)
d2bvals = eval.basis(1:365,daybasis365,Lfdobj=2)  # Lfdobj:linear differential operator 2阶线性微分算子

#循环建议不同number basis function 的基矩阵
for(i in c(1,2,3,6,9,12,15,20,26,52,104,182)){
  X = bvals[,1:(2*i+1)] # basis function构成的matrix
  yhat = X%*%solve(t(X)%*%X)%*%(t(X)%*%y)
  d2yhat = d2bvals[,1:(2*i+1)]%*%solve(t(X)%*%X)%*%(t(X)%*%y)  # 哪个函数？？
  
  
  # 绘制yhat拟合图像
  plot(1:365,y,col=2,cex=1.5,xlab='day',ylab='precipitation',
       cex.lab=1.5,cex.axis=1.5)
  lines(1:365,yhat,lwd=2,col=4)
  title(paste("No of basis", ncol(X)))
  
  # 绘制二阶导D2(X)图像，关注取值范围
  plot(1:365,d2yhat,col=4,type='l',lwd=2,cex.lab=1.5,cex.axis=1.5,
       ylab = 'D2 precipitation',xlab='day',main='Vancouver')
}

#### Simulation studies ####
## bias-variance simulation:
# (1) fit Vancouver precipitation by B-splines, to get x(ti) and pretend this is the `truth'. 
# (2) calculate `errors' ei = yi-x(ti), and create new `data' 
#     by randomly re-arranging the errors yi* = x(ti)+ei*
# (3) fit the new data using a Fourier basis 
#     and repeat nrep times to calculate bias and variance from sample. nrep = 10,100,1000.

# simulation重复次数
nrep=10
# 每月末最后一号作为knots, cumsum()累加
knots = cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31)) 
# B-spline method 建立basis function with norder = 6
bbasis = create.bspline.basis(rangeval=c(0,365),breaks=knots,norder=6)


tbvals = eval.basis(1:365,bbasis)
# 估计y
yhat = tbvals%*%solve( t(tbvals)%*%tbvals )%*%t(tbvals)%*%y
plot(x = 1:length(yhat), y = yhat, xlab = "Time in a Year /day", 
     ylab = "Observed precipitation values/mm", cex.lab=1.2,cex.axis=1.5)
abline(v=cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31)), col="gray",lwd=2)
special <- list(x=knots[-1], y = yhat[knots])
points(special, col="red", pch=16, cex=2)

# 对后续出现Runge现象跌落的预测进行重新估计
# yhat[351:365]存在震荡,变换之后的数据线性递减
# 不明确目的,可以不加
plot(yhat)
plot(yhat)[351:365]
yhat[351:365] = yhat[350] + (1:15)/16*(yhat[1]-yhat[350])

# 真实模型的error, 365个点估计
err = y-yhat


SimRes = array(0,c(365,nrep,50)) # 形成365*10矩阵, 50个 

for(i in 1:nrep){
  # 造新数据，给yhat从365个误差中随机加一个error
  ty = yhat + err[sample(365)] 
  
  for(j in 1:50){
    tbvals = bvals[,1:(2*j+1)]
    SimRes[,i,j] = tbvals%*%solve( t(tbvals)%*%tbvals )%*%t(tbvals)%*%ty 
    # 赋值：SimRes[,i,j]: 第i次模拟(第i列),(j block中)用2j+1个Fourier basis function得到的估计 
    #       SimRes[,i,1]: 表示使用3个Fourier basis function下, 每次重复实验计算得到的yhat new 存在每列中,
    #       也就是一个matrix (365*nrep) 是特定数目basis func条件下, 对yhat 重复nrep次的估计 
  }
}

## 计算var, bias, mse
Vars = apply(apply(SimRes,c(1,3),var),2,mean)  
# Simres三维数组, apply(SimRes,c(1,3),var)表示对同行(same day)同block(same basis func)求重复次数的方差
#                 dim 为365*50 (number of days)*(number of all possible basis func)

Bias = apply( (matrix(yhat,365,50)-apply(SimRes,c(1,3),mean))^2,2,mean)
# yhat 在simulation中表示真值, apply(SimRes,c(1,3),mean)表示给定i day,k basis条件下,average on repetition
#

MSE = apply( (array(yhat,c(365,nrep,50))-SimRes)^2,3,mean)
# SimRes 是\hat{X_ijk}, margin =3 说明对每个block I*J个element求均值


# 画图
par(mfrow=c(1,1))
stats = cbind(Vars,Bias,MSE)
matplot(1+2*(3:10),stats[3:10,],type='l',lwd=3,xlab='Number of Fourier Basis Functions',ylab='Estimation',main ="Bias-Variance Simulation, nrep = 10",cex.lab=1.2,cex.axis=1.5)
legend(16,0.080,c('Variance','Bias^2','MSE','Optimal number of Basis'),col=c(1,2,3,"blue"),lty=c(1,2,3,1),lwd=3)
abline(v=1+2*which.min(MSE),lty=2,lwd=1,col="blue")
text(13.25, 0.015, adj=0, 'min(MSE) with 13 Fourier Basis',col = "blue")


# install package:xtable to create table in latex
library(xtable)
library(stargazer)
stargazer(stats[1,10,])
xtable(stats[1:10,],digits=3,caption="Simulation Rusults with rnep = 10")


## minimize the cross-validated score:
smsse = rep(0,50)

for(i in 1:50){
  tbvals = bvals[,1:(2*i+1)]
  S = tbvals%*%solve( t(tbvals)%*%tbvals )%*%t(tbvals)  
  h = (1-diag(S))^2  # 为什么diag(S)各个元素都= 0.2767？
  errs = y - S%*%y   # yhat = S%*%y
  
  smsse[i] = sum( errs^2/h )  # sse of different number of basis functions
}  

# 绘制cv随基函数个数的图像(简化用SSE的表达)
plot(1+2*(2:20),smsse[2:20],type='b',lwd=3,xlab='number of Fourier Basis Functions',
     ylab='prediction error',
     main='Cross-Validated Error',cex.lab=1.5,cex.axis=1.5)
abline(v=1+2*which.min(smsse),lty=2,lwd=2,col="red")
text(22,377, adj=0, 'min(CV)=365.3: ')
text(22,375, adj=0, '21 Fourier Basis')
#abline(v=seq(5,40,5),col="gray",lwd=2)
#abline(h=seq(370,385,5),col="gray",lwd=2)
abline(v=13,lty=2,lwd=2,col="blue")
text(14,377, adj=0, 'min(MSE):')
text(14,375,adj=0," 13 Fourier Basis")
legend(35,370,c('min(SSE)','min(MSE)'),col=c("red","blue"),lty=c(2,2),lwd=3)


# 13 basis functions and look at standard errors:

tbvals = bvals[,1:13] # \phi(t)={1,sint,cost,...,sin6t,cos6t},# of basis = 13
S = tbvals%*%solve( t(tbvals)%*%tbvals )%*%t(tbvals)
yhat = S%*%y

err = (y-yhat)
sig = sum( (y-yhat)^2 )/(365-13) # sig = SSE/n-p,where df = n-p, p=13

Sig = sig*S%*%t(S) # variance matrix of X(t)
off = 2*sqrt(diag(Sig))

# 绘制拟合曲线及置信区间
plot(1:365,y,col=1,xlab='day in a year',ylab='precipitation/mm',cex.lab=1.5,cex.axis=1.5,
     main='Vancouver')
lines(1:365,yhat,lwd=2,col=2)
lines(1:365,yhat+off,col=4,lwd=2,lty=2)
lines(1:365,yhat-off,col=4,lwd=2,lty=2)
legend(35,8.5,c('estimation','95% CI'),col=c("red","blue"),lty=c(1,2),lwd=3)

