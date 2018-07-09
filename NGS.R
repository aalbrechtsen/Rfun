estPem<-function(like){ #numeric optimazition EM
    pk<-0.001 #start
    for(tal in 1:20){
      w0<-like[,1]*(1-pk)^2
      w1<-like[,2]*2*pk*(1-pk)
      w2<-like[,3]*pk^2
      pk<-mean((w1+2*w2)/(2*(w0+w1+w2)))
    }
    pk
  }


getLikes<-function(x,d=5,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}


emLog <- function(p=0.1,loglike,iter=30,accu=0.000001,accu2=0){
    numInds = ncol(loglike)
    temp_p = p
  for(it in 1:iter){
    sum=0
    for(i in 1:numInds){
      W0=exp(loglike[i*3+0-2])*(1-p)^2
      W1=exp(loglike[i*3+1-2])*2*p*(1-p)
      W2=exp(loglike[i*3+2-2])*p^2
      sum = sum + (W1+2*W2)/(2*(W0+W1+W2))
    }

    p=sum/numInds
    if((p-temp_p<accu & temp_p-p<accu) | (p/temp_p<1+accu2 & p/temp_p>1-accu2))
      break
    temp_p=p;
}
    p
}


estPem<-function(like){ #numeric optimazition EM
    pk<-0.001 #start
    for(tal in 1:20){
      w0<-like[,1]*(1-pk)^2
      w1<-like[,2]*2*pk*(1-pk)
      w2<-like[,3]*pk^2
      pk<-mean((w1+2*w2)/(2*(w0+w1+w2)))
    }
    pk
  }





if(FALSE){

ind<-1000
dep<-2
freqPop <- 0.1
genoTrue <- rbinom(ind,2,freqPop)
freqTrue <- mean(genoTrue)/2
l1<-getLikes(genoTrue,dep,norm=TRUE)

c(popFreq=freqPop,trueFreq=freqTrue,estFreq=estPem(l1),emLog=emLog(0.1,log(l1))




}
