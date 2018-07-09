freq.rr<-function (rr = 4.5, p = 0.15, kp = 0.1,freq=F){#kp prevalence
  f.mod <- c(rr^2, rr,1)
  q <- 1 - p
  fhw <- c(p^2, 2 * p * q, q^2)
  pi <- kp/sum(f.mod * fhw)
  if (pi <= 0 | pi >= 1) {
    warning("The combination of p, kp, and rr produces an unrealistic value of pi.")
    ret <- NA
  }
  else {
    fe <- rbind(fhw, fhw)
    dimnames(fe) <- list(c("Case", "Control"), c("AA", "Aa","aa"))
    f <- fe * rbind(f.mod * pi, 1 - f.mod * pi)
    ret<-f/ apply(f, 1, sum)
    if(freq)
      ret<-(ret[,2]+2*ret[,1])/2
  }
  ret
}
