makelookuparray <- function(CA,S,coeff,B,Ls,N,J,a,Hj) {
  #create lookup table for jam with lower gap
  Hpts=seq(0,4,0.01) #HJ/HBF
  rivarray <- array(0,dim=c(length(Hpts),6,N))
  Npts=seq(1,N,1)
  Cp0=2/3
  Cb=Cp0*coeff/S-1
  for(i in Npts){
    colhj <- Hpts
    colQbelow <- (colhj<a[i])*B[i]*sqrt(S[i]/coeff[i]*9.8*colhj^3)
    colQjam <- (colhj<=Hj[i])*(colhj>=a[i])*B[i]*((2*9.8*(colhj-a[i])^3/(3^(3/2)*CA[i]))^(1/2)+(Cp0/(1+Cb[i]*a[i]/colhj)*9.8*a[i]^2*colhj)^(1/2))
    colQjam[is.nan(colQjam)] = 0 
    colQabove <- (colhj>Hj[i])*(2*sqrt(2*9.8)/3*B[i]*(colhj-Hj[i])^(3/2)+max(colQjam)) #value when jam is full
    colQabove[is.nan(colQabove)] = 0
    colQ <- colQbelow+colQjam+colQabove
    colh0 <- ((colQ/B[i])^2*coeff[i]/(9.8*S[i]))^(1/3)
    colV <- Ls[i]*B[i]*colh0+B[i]/(2*S[i])*(colh0-colhj)^2*J[i]
    colA <- B[i]*colh0+B[i]/(2*S[i]*Ls[i])*(colh0-colhj)^2*J[i]
    colhe<- colh0+1/(2*S[i]*Ls[i])*(colh0-colhj)^2*J[i] #effective constant h should be able to replace
    rivarray[,,i] = c(colQ,colh0,colhj,colV,colA,colhe)
}
  return(rivarray)
}
  