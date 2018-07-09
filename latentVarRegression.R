
fit_lm <- function(design, pheno, weights=NULL){
  if (dim(design)[1] != length(pheno)) { stop("dimention mismatch in fit_lm") }
  if (is.null(weights)){
    weights <- rep(1, length(pheno))
  }
  z <- (weights == 0)
  xx <- design[!z, , drop=FALSE]
  yy <- pheno[!z]
  ww <- weights[!z]
  ttol <- 1e-07
  wwts <- sqrt(ww)
  #z <- .Call(stats:::C_Cdqrls, xx * wwts, yy * wwts, ttol, FALSE)
  xw <- xx * wwts
  yw<- yy * wwts
  coef <- (solve(t(xw) %*% xw)) %*% (t(xw) %*% yw)
  df <- length(yy) - length(coef)
  residuals <- pheno - design %*% coef
  return(list(coefficients = coef, df.residual = df, residuals = residuals))
}


#' Calculate joint probabilities of phenotypes and ancestral states for quantitative traits
#'
#' @param params Parameters, first genetic effects, then intercept and additional covariate effects and last sd (length 2/3 + p + 1)
#' @param phenos Individual phenotypes (augmented to length 4n, i.e. each repeated 4 times)
#' @param design Augmented design matrix ( matrix 4n x length(params)-1 )
#' @param prior Probability of ancestral state conditional on genotype (vector length 4n) - see ?make_mixture_props
#' @return Joint probability of phenotype and state given parameters (matrix 4 x n)
make_joint_lm <- function(params, phenos, design, prior) {
  probs<-dnorm(x = phenos,
               mean = design%*%(params[-length(params)]),
               sd = params[length(params)])
  matrix(probs * prior,nr=4)
}

#' Calculate log of joint probabilities of phenotypes and ancestral states for quantitative traits
#'
#' @param params Parameters (first genetic effects, then intercept and additional covariate effects and last sd)
#' @param phenos Individual phenotypes augmented to length 4 x N
#' @param design Augmented design matrix, 4N rows
#' @return Matrix with 4 rows and N columns, with log of joint probability of phenotype and state given parameters
make_joint_lm_log <- function(pars, phe, desi, prior,nLatent=3) {
  probs <- dnorm(x = phe,
                 mean = desi%*%(pars[-length(pars)]),
                 sd = pars[length(pars)], log=T)
  matrix(probs + log(prior), nr=nLatent)
}





#' colSums in log-space
colSumsLog <- function(x){
  k <- apply(x, 2, max)
  log(rowSums(exp(t(x) - k))) + k
}

#' EM update step for quantitative traits, calculations done in log space
#'
#' @param pars Model parameters (length 2/3 + p +1 )
#' @param phen Individual phenotypes augmented to length 4 x N
#' @param dsgn Augmented design matrix (matrix 4n x length(params)-1)
#' @param prior Probability of ancestral state conditional on genotype (vector length 4n) - see ?make_mixture_props
#' @return list of two: pars_new are new parameter estimates (effects and sd) and mll_new is the updated minus log likelihood
update_em_lm_log<-function(pars, phen, dsgn, prior,nLatent=3) {
  ## joint pheno and state posterior prob
  joint_pheno_state_log <- make_joint_lm_log(pars = pars, phe = phen, desi = dsgn, prior = prior)
  prob_pheno_log <- colSumsLog(joint_pheno_state_log)
  mll_new<- -sum(prob_pheno_log)
  df <- length(phen)/nLatent-ncol(dsgn)
  wgts_log<- c(t(t(joint_pheno_state_log) - prob_pheno_log))
  wgts <- exp(wgts_log)
  ## New parameter guess
  fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
  #fit <- lm(phen~.-1,weights = wgts,data = data.frame(phen,dsgn))
  sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
  pars_new <- c(fit$coefficients, sigm)

  fit <- fit_lm(design = dsgn, pheno = phen, weights = wgts)
  sigm <- sqrt(sum((fit$residuals^2) * wgts) / df)
  pars_new <- c(fit$coefficients, sigm)

  return(list(pars_new = pars_new, mll_new = mll_new))
}


#' EM update step wrapper for choosing quantitative or dichotomuous traits
update_em <- function(pars, phen, dsgn, prior, quant) {
  if (quant == TRUE) {
    #update_em_lm(pars, phen, dsgn, prior)
    update_em_lm_log(pars, phen, dsgn, prior)
  } else {
    update_em_bin(pars, phen, dsgn, prior)
  }
}


#' EM update step for dichotomuous traits
#'
#' @return list of two: pars_new are new parameter estimates (effects) and mll_new is the updated minus log likelihood
update_em_bin <- function(pars, phen, dsgn, prior) {
  ## joint pheno and state posterior prob
  joint_pheno_state <- make_joint_bin(pars = pars, phe = phen, desi = dsgn, prior = prior)
  ## Marginal probs of y (state summed out)
  prob_pheno <- colSums(joint_pheno_state)
  ## minus log likelihood (for old parameters actually)
  mll_new <- -sum(log(prob_pheno))
  ## Weights are the posterior state probabilities given pheno and pars
  wgts <- c(t(t(joint_pheno_state) / prob_pheno))
  ## New parameter guess
  fit <- fit_bin(design = dsgn, pheno = phen, weights = wgts)
  pars_new <- fit$coefficients
  return(list(pars_new = pars_new, mll_new = mll_new))
}


#' EM algorithm controller
#'
#' @param initial Start values for optimization (if quantitative length is ncol(desi)+1 else length is ncol(desi))
#' @param maxI Max number iterations of EM algo
#' @param phe Observed phenotypes 4n long
#' @param desi Sesign matrix 4n times 2+nCov
#' @param pri Prior dist over states, 4n long
#' @param qua Is trait quantitative? true/false
#' @param tole Convergence tolerence
#' @return list of par (estimates), value (minus log likelihood), counts (iterations), convergence and about (how did algo terminate)


control_em <- function(phe, desi, pri, qua, tole=1e-4,start, maxI=100){

    if(missing(start)){
        start <- rnorm(ncol(desi))
        if(qua)
            start<-c(start,sd(phe))
    }

    
    pars_old <- start
    mll_old <- Inf
    for (iter in 1:maxI) {
        update <- update_em(pars = pars_old, phen = phe, dsgn = desi, prior = pri, quant=qua)
        pars_new <- update$pars_new
        mll_new <- update$mll_new
        if (mll_new > mll_old + 10e-6) {
            print(mll_new)
            print(mll_old)
            em_out <- list(par=pars_new, value = mll_new, counts = c(iter, iter), convergence = 1, about = "mll_new > mll_old + 10e-6")
            print("EM step in wrong direction.")
            break
        } else if (mll_old - mll_new < tole) {
            em_out <- list(par = pars_new, value = mll_new, counts = c(iter, iter), convergence = 0, about = "mll_old - mll_new < tol")
            break
        } else if ( iter == maxI){
            em_out <- list(par = pars_new, value = mll_new, counts = c(iter, iter), convergence = 1, about="iter == maxI")
        } else {
            pars_old <- pars_new
            mll_old <- mll_new
        }
    }
  return(em_out)
}




lmLatent<-function(y,xp,z){

    if(missing(z)){
       # l<-lm(y~1)
        Xe<-cbind(rep(1,length(y)))
    }
    else{
        z<-cbind(z)
      #  l<-lm(y~1+z)
        Xe<-cbind(1,z) 
    }
    if(is.null(dim(xp))){
        mat<-matrix(0,nrow=3,ncol=length(y))
        mat[xp+(1:length(y)-1)*3+1]<-1
        xp<-t(mat)
    }
    
    N<-length(y)
    nLatent<-nrow(xp)
    phe <- rep(y,each=nLatent)
    pri <- as.vector(xp)
    qua <- TRUE
    nPara <- ncol(Xe)
    desi <- rep(Xe,each=nLatent)
    dim(desi) <- c(nLatent*N,nPara)
    
    desi <- cbind(rep(1:nLatent-1,N),desi)
    
  control_em(phe, desi, pri, qua)
  
}




normScoreCov<-function(y,xp,z){
  if(missing(z)){
    l<-lm(y~1)
    Xe<-cbind(rep(1,length(y)))
  }
  else{
    z<-cbind(z)
    l<-lm(y~1+z)
    Xe<-cbind(1,z) 
  }
    if(is.null(dim(xp))){
      mat<-matrix(0,nrow=3,ncol=length(y))
      mat[xp+(1:length(y)-1)*3+1]<-1
      xp<-t(mat)
    }
    else if(nrow(xp)==3)
      xp<-t(xp)
  Ex<-as.vector(xp%*%0:2)
  Ex2<-as.vector(xp%*%(0:2)^2)
#  yfit<-X%*%MASS::ginv(t(X)%*%X)%*%t(X)%*%y
  yTilde<-l$fitted.values
  var<-sum(l$residuals^2)/(length(y)-2)

  Vb<-0
  Vaa<-matrix(0,ncol=ncol(Xe),nrow=ncol(Xe))
  Vab<-matrix(0,ncol=1,nrow=ncol(Xe))
  Vbb<-0
  U<-0
  for(tal in 1:length(y)){
   U<-U+(y[tal]-yTilde[tal])/var*Ex[tal]
   Vaa<-Vaa+1/var*Xe[tal,]%*%t(Xe[tal,])
   Vab<-Vab+1/var*Ex[tal]*cbind(Xe[tal,])
   k<-(y[tal]-yTilde[tal])^2/var^2
   Vbb<-Vbb+(1/var-k)*Ex2[tal]+k*Ex[tal]^2


  }
  I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  1-pchisq(U^2/I,1)
}



normScore<-function(y,xp){
    xp <- t(xp)
  if(is.null(dim(xp))){
    mat<-matrix(0,nrow=3,ncol=length(y))
    mat[xp+(1:length(y)-1)*3+1]<-1
    xp<-t(mat)
    }
    
  Ex<-as.vector(xp%*%0:2)
  Ex2<-as.vector(xp%*%(0:2)^2)
  yTilde<-mean(y)
  var<-var(y)
  
  U<-sum((y-yTilde)/var*Ex)
  Vaa<-1/var*length(y)
  Vab<-1/var*sum(Ex)
  k<-(y-yTilde)^2/var^2
  Vbb<-sum( (1/var-k)*Ex2+k*Ex^2)
  I<-Vbb-Vab^2/Vaa
  1-pchisq(U^2/I,1)
}
binScoreEnv<-function(y,xp,z=NULL){
  #not done
    if(is.null(dim(xp))){
    mat<-matrix(0,nrow=3,ncol=length(y))
    mat[xp+(1:length(y)-1)*3+1]<-1
    xp<-t(mat)
    }
    if(is.null(z)){
      A<-cbind(rep(1,length(y)))
      d<-1
    }
    else{
      A<-cbind(rep(1,length(y)),z)
      d<-dim(A)[2]
    }

    var<-1
    Ex<-as.vector(xp%*%0:2)
    Ex2<-as.vector(xp%*%(0:2)^2)
    yTilde<-A%*%solve(t(A)%*%A)%*%t(A)%*%y
    #yTilde<-mean(y)
    U<-sum((y-yTilde)/var*Ex)

    Vaa<-0
    Vba<-0
    Vbb<-0
    for(i in 1:length(y)){
      Vaa<-Vaa+yTilde[i]*(1-yTilde[i])*A[i,]%*%t(A[i,])
      Vba<-Vba+yTilde[i]*(1-yTilde[i])*A[i,]*Ex[i]

    }
    Vba<-t(t(Vba))
    Vbb<-sum((yTilde*(1-yTilde)-(y-yTilde)^2)*Ex2+(y-yTilde)^2*Ex^2)

    I<-Vbb-t(Vba)%*%solve(Vaa)%*%Vba
    1-pchisq(U^2/I,1)
}

binScore<-function(y,xp){
    xp <- t(xp)
   if(is.null(dim(xp))){
    mat<-matrix(0,nrow=3,ncol=length(y))
    mat[xp+(1:length(y)-1)*3+1]<-1
    xp<-t(mat)
    }

    Ex<-as.vector(xp%*%0:2)
    Ex2<-as.vector(xp%*%(0:2)^2)
    yTilde<-mean(y)
 
    U<-sum((y-yTilde)*Ex)
    Vaa<-yTilde*(1-yTilde)*length(y)
    Vab<-yTilde*(1-yTilde)*sum(Ex)
    Vbb<-sum((yTilde*(1-yTilde)-(y-yTilde)^2)*Ex2+(y-yTilde)^2*Ex^2)
    I<-Vbb-Vab^2/Vaa
    1-pchisq(U^2/I,1)
}



lmFull<-function(y,xp){
  y<-y
  if(is.null(dim(xp))){
    mat<-matrix(0,nrow=3,ncol=length(y))
    mat[xp+(1:length(y)-1)*3+1]<-1
    xp<-mat
  }

  
 
  lmLike<-function(x){
    sd<-x[1]
    alpha<-x[2]
    beta<-x[3]
    l0<-dnorm(y,sd=sd,mean=alpha)
    l1<-dnorm(y,sd=sd,mean=alpha+beta)
    l2<-dnorm(y,sd=sd,mean=alpha+2*beta)
    L<--sum(log(xp[1,]*l0+xp[2,]*l1+xp[3,]*l2))
    L
  }
  lmLike0<-function(){
    sd<-sd(y)
    alpha<-mean(y)
    L<--sum(log(dnorm(y,sd=sd,mean=alpha)))
    L
  }

  par<-optim(c(sd(y),rnorm(1),rnorm(1)),lmLike,method="L-BFGS-B",lower=c(sd(y)/5,-1,-1),upper=c(sd(y),2,2))
  null<-lmLike0()
  names(par$par) <- c("sigma","intersept","beta")
  par$pval <- 1-pchisq(2*(null-par$value),1)
  par
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



getLikes<-function(x,d=5,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}
##x is a N*M matrix
##y is a M vector
## each column is divided by a entry in M
matScale <- function(x,y,type=1){
    if(type==1)
        x %*% diag(1/y)
    else if(type==2)
        x/y[col(x)]

}

getPost <- function(GL,f){
    PP <- GL * c((1-f)^2,2*f*(1-f),f^2)
    matScale(PP,colSums(PP),type=2)
}


if(FALSE){
    source("latentVarRegression.R")

    N <- 100000
    f <- 0.2
    sd=0.4
    G<-rbinom(N,2,f)
    beta <- 0.75
    beta <- 0.015
    ## beta <- 0
    y<-rnorm(N,sd=sd) + G*beta


    ## 3 x N
    GL <- getLikes(G,d=3,norm=T)
    ## GL <- getLikes(G,d<-1+5*rank(y)/length(y),norm=T)


    ## est freq
    f.ml<- estPem(GL)#EM

    ## post
    PP <- getPost(GL,f)

    ## dosage
    EG <- PP[2,]+2*PP[3,]
    
    ## best
    best<-(PP[1,]<PP[2,]|PP[1,]<PP[3,])+(PP[2,]<PP[3,]&PP[1,]<PP[3,])

    #TRUE
    coef(lg <- lm(y~G))["G"]
    #dosage
    coef(le <- lm(y~EG))["EG"]
    #Weighted reg
    coef(lw <- lm(rep(y,each=3)~rep(0:2,N),weights=as.vector(PP)))[2]
    ## highest GL (calling)
    coef(lb <- lm(y~best))["best"]
    ## full latent
    system.time(    ll <- lmFull(y,PP));    ll$par["beta"] 
    ## full latent EM
    system.time(ll2 <- lmLatent(y,PP))

    getPval <- function(x) summary(x)$coefficients[2,4]









  
}
