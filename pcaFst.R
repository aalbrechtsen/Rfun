## PCA (Eigenstrat/smartpca style) and Weir & Cockerham (1984) Fst for genotype matrices.
## Moved from the popgen.dk software wiki ("Rscripts"). See pcaFst.md for examples.

eigenstrat<-function(geno,maxMis=0,minMaf=0.01){
## geno: snp x ind matrix of genotypes \in 0,1,2
##maxMis maximum allowed missing genotypes for a site

  nMis <- rowSums(is.na(geno))
  freq <- rowMeans(geno,na.rm=T)/2               # get allele frequency
  keep <- freq>minMaf&freq<(1-minMaf) & nMis<=maxMis         # remove sites with non-polymorphic data
  freq<-freq[keep]
  geno<-geno[keep,]
  snp<-nrow(geno)                           #number of snps used in analysis
  ind<-ncol(geno)                           #number of individuals used in analuysis
  M <- (geno-freq*2)/sqrt(freq*(1-freq))       #normalize the genotype matrix
  M[is.na(M)]<-0
  X<-t(M)%*%M                               #get the (almost) covariance matrix
  X<-X/(sum(diag(X))/(snp-1))
  E<-eigen(X)

  mu<-(sqrt(snp-1)+sqrt(ind))^2/snp         #for testing significance (assuming no LD!)
  sigma<-(sqrt(snp-1)+sqrt(ind))/snp*(1/sqrt(snp-1)+1/sqrt(ind))^(1/3)
  E$TW<-(E$values[1]*ind/sum(E$values)-mu)/sigma
  E$mu<-mu
  E$sigma<-sigma
  E$nSNP <- nrow(geno)
  E$nInd <- ncol(geno)
  class(E)<-"eigenstrat"
  E
}

plot.eigenstrat<-function(x,col=1,...)
  plot(x$vectors[,1:2],col=col,...)

print.eigenstrat<-function(x,...)
  cat("statistic",x$TW,"\n")


## x: ind x snp matrix of genotypes (0,1,2); pop: population label for each individual
WC84<-function(x,pop){
  x<-x[,!apply(is.na(x),2,any)]

  popNam<-sort(unique(pop))

  #number ind each population
  n<-table(pop)
  #number of populations
  npop<-nrow(n)
  #average sample size of each population
  n_avg<-mean(n)
  #total number of samples
  N<-length(pop)
  #frequency in samples
  p<-t(sapply(popNam,function(y) colMeans(x[pop==y,])/2))
  #average frequency in all samples (apply(x,2,mean)/2)
  p_avg<-as.vector(n%*%p/N )
  #the sample variance of allele 1 over populations
  s2<-1/(npop-1)*(apply(p,1,function(x){((x-p_avg)^2)})%*%n)/n_avg
  #average heterozygote frequency for allele 1
  h_avg<-apply(x==1,2,sum)/N
  #nc (see page 1360 in wier and cockerhamm, 1984)
  n_c<-1/(npop-1)*(N-sum(n^2)/N)
  #variance betwen populations
  a <- n_avg/n_c*(s2-(p_avg*(1-p_avg)-(npop-1)*s2/npop-h_avg/4)/(n_avg-1))
  #variance between individuals within populations
  b <- n_avg/(n_avg-1)*(p_avg*(1-p_avg)-(npop-1)*s2/npop-(2*n_avg-1)*h_avg/(4*n_avg))
  #variance within individuals
  c <- h_avg/2

  #inbreedning (F_it)
  F <- 1-c/(a+b+c)
  #(F_st)
  theta <- a/(a+b+c)
  #(F_is)
  f <- 1-c/(b+c)
  #weigted average of theta
  theta_w<-sum(a)/sum(a+b+c)
  list(F=F,theta=theta,f=f,theta_w=theta_w,a=a,b=b,c=c,total=c+b+a)
}
