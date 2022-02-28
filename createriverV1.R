#set up river network to be called by findoutflowV1.R
createriverV1 <- function() {

g=9.8

#Properties from Usway Burn @ Shillmoor (Hey & Thorne Site 57)
S=0.008479
HBF=0.78 #m
BBF=9.1 #m
Ds50=0.1135 #m
a=0.5*HBF #gap height m
HJ=HBF #jam height m

Cf0=(5.75*log10(2*HBF/Ds50))^-2  #Julien 1998
QBF=BBF*sqrt(g*HBF^3*S/Cf0) 

maxQ=QBF
minQ=0.2*maxQ

#define high, med low density barrier aCA=0.8, 0.5, 0.2
aCA=0.2 
CA=2*Cf0/(3^(3/2)*aCA^2*S)
HJ_QBF=sqrt(3)*(CA*(QBF/BBF)^2/(2*g))^(1/3)
LBW_QBF=(HJ_QBF-HBF)/S

#spacing density: LBW_BF=100%, 50%, 20%, 10% of Ls b=1,2,5,10
b=1
Ls=b*LBW_QBF

N=100; #Number of segments

#river parameters  
LR = Ls*(N+1) #because Nsegs=Lx-2 in model--wanted reach length of computation to be 15000 
Lseg = Ls
Nsegs = N
slope=S

#for user to set N=100,50,25,10 jams
rep100=c(1)
repN100=N
rep50=c(0,1)
rep50=c(1,0)
repN50=N/2
rep25=c(0,0,0,1)
rep25=c(1,0,0,0)
repN25=N/4
rep10=c(0,0,0,0,0,0,0,0,0,1)
rep10=c(1,0,0,0,0,0,0,0,0,0)
repN10=N/10
riv <- list(CA=c(0,rep(CA,Nsegs),0))
riv["J"] <-list(J=c(0,rep(rep25,repN25),0)) 
#J and CA are slightly redundant; of course if CA=0 then it doesn't matter what J is and vice versa. 
#But, it seems a bit simpler for the user to explicitly state that a jam is not present (J=0), rather than have to set CA=0

riv["x"] <- list(x=c(0,seq(Lseg,(Nsegs*Lseg),Lseg),LR))
riv["z"] <- list(z=slope*max(riv$x)-slope*riv$x+0)
riv["H"] <- list(H=rep(HBF,length(riv$x)))
riv["B"] <- list(B=rep(BBF,length(riv$x)))
riv["S"] <- list(slope=rep(S,length(riv$x)))
riv["QBF"] <- QBF
riv["coeff"] <- list(coeff=rep(Cf0,length(riv$x)))
riv["a"] <- list(a=rep(a,length(riv$x)))
riv["HJ"] <- list(HJ=rep(HJ,length(riv$x)))

return(riv)
}
