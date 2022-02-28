#sources files and runs findoutflowV1()
#enter path below to source files
pathstr <- "" 
source(paste(pathstr,"getinflowV1.R",sep=""))
source(paste(pathstr,"createriverV1.R",sep=""))
source(paste(pathstr,"makelookuparray.R",sep=""))
source(paste(pathstr,"odes3.R",sep=""))
source(paste(pathstr,"findoutflowV1.R",sep=""))

findoutflowV1()
