#system of ODEs comprising network model
odes3 <- function(t,y,parms) {
  B <- parms[[1]]
  S <- parms[[2]]
  coeff <- parms[[3]]
  g <- parms[[4]]
  CA <- parms[[5]]
  J <- parms[[6]]
  tInR <- parms[[7]]
  QInR <- parms[[8]]
  N <- parms[[9]]
  ad <- parms[[10]]
  l <- parms[[11]]
  a <- parms[[12]]
  Hj <- parms[[13]] #height of top of jam
  rivarray <- parms[[14]]
  
  A=y;
  Npts <- seq(1,N,1)
  harray  <-array(0,c(1,N))
  Qarray  <-array(0,c(1,N))  
  for(i in Npts){
    colchoice<-J[i]+2
    harray[i]=interp1(rivarray[,5,i],rivarray[,colchoice,i],A[i]) #A,he or h0, hj
    Qarray[i]=interp1(rivarray[,5,i],rivarray[,1,i],A[i]) #A,Q
  }
  
  h=t(harray)
  Q=t(Qarray)
  Qin=interp1(tInR*3600, QInR, t)
  q=c(1,rep(0,N-1))*Qin
  dydt=l^(-1)*(t(ad) %*% Q - Q + q)  
  return(list(as.vector(dydt)))
}  
