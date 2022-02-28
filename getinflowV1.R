getinflowV1 <- function(QBF){
#set inflow to first segment of river channel
maxQ = QBF;
minQ = 0.2*maxQ;
mu = 6;
sigma = 1;
xi=seq(0,24*10,0.05) #hr
yi=1/(sigma*sqrt(2*pi))*exp(-0.5*((xi-mu)/sigma)^2)
coeff=(maxQ-minQ)/max(yi);
yi=yi*coeff+minQ; 
inflow <- list(t=xi) 
inflow["Q"] <- list(Q=yi) #m3/s
return(inflow)
}

