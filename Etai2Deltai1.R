Etai2Deltai1 <- function(etai, gamma, lambda, sigma, q.points, Z, W){ 
  ## Compute delta for shifted gamma model
  expit <- function(aa){exp(aa)/(1+exp(aa))}
  n     <- length(etai)
  #print(etai)
  gam   <- gamma
  if (length(gam)==1){ gam<-rep(gamma, n) }
  lam   <- lambda
  if (length(lam)==1){ lam<-rep(lambda, n) }
  sig   <- sigma
  if (length(sig)==1){ sig<-rep(sigma, n) }
  
  if (length(sig) != n | length(gam) !=n) { stop("Error in etai.2.deltai") }

  deltai  <- rep(NA, n)
  deltai  <-.Call("DeconvolveGL_CALL", etai, gam, lam, sig, Z, W)
  deltai
 # print(deltai)
}
