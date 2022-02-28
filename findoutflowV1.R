#computes and plots outflow from a 1D river network with series of jams
#input channel, inflow parameters in createriverV1, getinflowV1
#calls createriverV1.R, getinflowV1.R, makelookuparray.R, odes3.R
#plots inflow & outflow hydrographs and max change in level at first section

findoutflowV1 <- function(coeff=1,jams=1){
#load necessary libraries
library(pracma)  
library(deSolve)
library(Matrix)
  
#import river parameters specified in createriverV1() and inflow q,t from getinflowV1()
#jams options allows running the same network with and without a series of jams
#input coeff allows input of bed resistance, currently C_f
if (jams==1){
  riv <- createriverV1()
  inflow <-getinflowV1(riv$QBF)
} else{
  riv <- createriverV1()
  inflow <-getinflowV1(riv$QBF)
  riv$J <- 0*riv$J #same network as jams but with J=0 at all nodes
  riv$coeff <- coeff*(1/riv$coeff[1])*riv$coeff #gives array with coeff=input in function call
}

#add 1h to inflow time--avoids apparent bug where t queried by lsodes > max(inflow$t) 
inflow$t <- append(inflow$t,inflow$t[length(inflow$t)]+1)
inflow$Q <- append(inflow$Q,inflow$Q[length(inflow$Q)])

#set parameters for use in solving system
lx <- length(riv$x) #used in later lines
N <- lx-2; #number of segments--follows code accompanying Hankin et al. (2020)

ad <- sparseMatrix(i=1:(N-1),j=2:N,x=1,dims=c(N,N)) #set up for 1 single channel; would need to alter for sub-channels in network
l <- matrix(riv$x[2:(lx-1)]-riv$x[1:(lx-2)],N,1)
S <- matrix(-(riv$z[2:(lx-1)]-riv$z[1:(lx-2)])/(riv$x[2:(lx-1)]-riv$x[1:(lx-2)]),N,1)
CA <- matrix(riv$CA[2:(lx-1)],N,1) 
g <- matrix(rep(9.8,N),N,1)
B <- matrix(riv$B[2:(lx-1)],N,1)
coeff <- matrix(riv$coeff[2:(lx-1)])
J <- matrix(riv$J[2:(lx-1)],N,1)
a <- matrix(riv$a[2:(lx-1)],N,1)
HJ <- matrix(riv$HJ[2:(lx-1)],N,1)
rivarray <-makelookuparray(CA,S,coeff,B,l,N,J,a,HJ)  

#solve system
#1. constant inflow: find initial conditions (ICs)
parmsIC <- list(B,S,coeff,g,CA,J,inflow$t,rep(inflow$Q[1],length(inflow$t)),N,ad,l,a,HJ,rivarray)
baseflowguess=(B[1]*(riv$coeff[1]*(min(inflow$Q)/B[1])^2/(g[1]*S[1]))^(1/3))
yIC=rep(baseflowguess,N) #initial guess--could improve based on uniform flow
outIC <- lsodes(y=yIC,times=c(0.25,48*3600),func=odes3,parms=parmsIC)
#2. variable inflow: use ICs from (1)
parms<-list(B,S,coeff,g,CA,J,inflow$t,inflow$Q,N,ad,l,a,HJ,rivarray)
times=seq(0, (9*24)*3600,0.25*3600) #related to bug mentioned in Line 12-input only original time which reduces maximum t input by lsodes
outA <- lsodes(y=as.vector(outIC[2,2:(lx-1)]),times=times,func=odes3,parms=parms)

#find h,Q only for last node, N, to reduce time
#save in allout structure
allout <- list(A = outA[,N])
lenTimes=length(times)
hOut=rep(0,lenTimes)
QOut=rep(0,lenTimes)

for(i in seq(1,lenTimes)){
  colchoice<-J[N]+2
  disp(colchoice)
  hOut[i]=interp1(rivarray[,5,N],rivarray[,colchoice,N],outA[i,N])
  QOut[i]=interp1(rivarray[,5,N],rivarray[,1,N],outA[i,N])
}
allout["h"] <- list(h = hOut)
allout["Q"] <- list(Q = QOut)
allout["t"] <- list(t = times)
allout["Q_in"] <- list(Q_in = inflow$Q)
allout["t_in"] <- list(t_in = inflow$t*3600)

#calculations for plots
tpeakIn=allout$t_in[which.max(allout$Q_in)]/3600
tpeakOut=allout$t[which.max(allout$Q)]/3600
hIn=interp1(rivarray[,1,1],rivarray[,2,1],allout$Q_in)
hInJ=interp1(rivarray[,1,1],rivarray[,3,1],allout$Q_in)
CAinds=which(riv$CA != 0)

#plot results: input and output Q,t &
#max change in level at first segment
plot(allout$t_in/3600,allout$Q_in,type='l',col='light blue',xlab="time (h)",
     ylab=expression(paste("Discharge ","(m"^"3","/s)")),ylim=c(0,max(allout$Q_in)),xlim=c(0,1*24),
     main=paste("peak reduction",signif(max(allout$Q_in)-max(allout$Q),digits=2),"(m³/s)\n increase in time to peak",signif(tpeakOut-tpeakIn,4),'hr'))
lines(allout$t/3600,allout$Q,col='blue')
legend(15,5,legend=c("inflow","outflow"),col=c("light blue","blue"), lty=c(1,1), cex=0.8)
plot(allout$t_in/3600,hInJ/riv$H[1],type='l',col='brown',xlab='time (h)',
     ylab='water depth / bankfull depth', 
     main=paste('max change in level at first section',signif(max(hInJ-hIn),2),'m \n barrier C_A =',signif(riv$CA[CAinds[1]],3)),ylim=c(0,max(hInJ/riv$H[1])),xlim=c(0,1*24))
lines(allout$t_in/3600,hIn/riv$H[1],col="deepskyblue")
legend(10,1,legend=c("barrier","no barrier"),col=c("brown","deepskyblue"), lty=c(1,1), cex=0.8)
return(allout)
}
