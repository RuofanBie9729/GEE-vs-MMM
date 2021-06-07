GEE.var.kc3 <- function(data, model, betahat, phi=1, alpha){
### function to compute KC correction sd for GEE

data <- data[order(data$id),]
  ##
  id <- data$id
  n <- nrow(data)
  id <- data$id
  ids <- unique(id)
  K   <- length(ids) # number of clusters
  beta.hat <- betahat
  p   <- length(beta.hat)
  K_frac <- K/(K-p)
  X <- model.matrix(model, data)
  Y <- data$binY
  Y_k <- by(Y, data$id, function(y){as.vector(y)})
  X_k <- by(X, data$id, function(y){as.matrix(y)})
  eta_k <- lapply(1:K, function(x){X_k[[x]]%*%beta.hat})
  mu_k <- lapply(1:K, function(x){as.vector(expit(eta_k[[x]]))})
  Y_mu_k <- lapply(1:K, function(x){Y_k[[x]] - mu_k[[x]]})
  n_k <- lapply(1:K, function(x){length(mu_k[[x]])})
  A_k <- lapply(1:K, function(x){diag(n_k[[x]])*(mu_k[[x]]*(1-mu_k[[x]]))})
  D_k <- lapply(1:K, function(x){A_k[[x]]%*%X_k[[x]]})
  D_k_t <- lapply(D_k, "t") 
  S_k <- lapply(1:K, function(x){(phi)^(1/2)*diag(n_k[[x]])*sqrt((mu_k[[x]]*(1-mu_k[[x]])))})
  Sk_inv <- lapply(1:K, function(x){(phi)^(-1/2)*diag(n_k[[x]])*(mu_k[[x]]*(1-mu_k[[x]]))^(-1/2)})
  alph_mat_init <- lapply(1:K, function(x){matrix(alpha, nrow=n_k[[x]], ncol=n_k[[x]])-diag(n_k[[x]])*alpha})
  alph_mat <- lapply(1:K, function(x){alph_mat_init[[x]]+diag(n_k[[x]])})
  Vk <- lapply(1:K, function(x){S_k[[x]]%*%alph_mat[[x]]%*%S_k[[x]]})
  Vk_inv <- lapply(1:K, function(x){Sk_inv[[x]]%*%solve(alph_mat[[x]])%*%Sk_inv[[x]]})
  Wts_k  <- by(data$wts, data$id, function(y){y}) 
  D <- Reduce("rbind", D_k)
  tD <- t(D)
  V_inv <- as.matrix(bdiag(Vk_inv))
  Y <- as.vector(Y)
  mu <- as.vector(Reduce(c, mu_k))
  Y_mu <- diag(n)*(Y - mu)
  ##
  U_t <- tD%*%V_inv%*%Y_mu
  W2K <- lapply(1:K, function(x){matrix(Wts_k[[x]][1], nrow=n_k[[x]], ncol=n_k[[x]])})
  W2 <- as.matrix(bdiag(W2K))
  W <- diag(n)*as.vector(Reduce(c, Wts_k))
  pi <- lapply(1:K, function(x){1/Wts_k[[x]][1]})
  pi_sq <- lapply(1:K, function(x){pi[[x]]^2})
  delta <- lapply(1:K, function(x){(pi[[x]] - pi_sq[[x]])/pi[[x]]})
  Delta <- lapply(1:K, function(x){matrix(delta[[x]], nrow = n_k[[x]], ncol = n_k[[x]])})
  P <- lapply(1:K, function(x){matrix(pi[[x]], nrow = n_k[[x]], ncol = n_k[[x]])})
  P_sq <- lapply(1:K, function(x){matrix(pi_sq[[x]], nrow = n_k[[x]], ncol = n_k[[x]])})
  ##
  prob_k <- lapply(1:K, function(x){1/Wts_k[[x]]})
  prob <- as.vector(Reduce(c, prob_k))
  Pij_1 <- prob %*% t(prob)
  ## this should be a block diag matrix
  Pij_2 <- as.matrix(bdiag(P_sq))
  Pij_3 <- Pij_1 - Pij_2 
  Pij_4 <- as.matrix(bdiag(P))
  Pij <- Pij_3 + Pij_4 
  Delta_tild <- as.matrix(bdiag(Delta))
  Delta_tild_e <- Delta_tild-as.matrix(bdiag(Delta))
  Delta_tild_c <- as.matrix(bdiag(Delta))
  H <- Reduce("+", lapply(1:K, function(x){Wts_k[[x]][1]*D_k_t[[x]]%*%Vk_inv[[x]]%*%D_k[[x]]}))
  H_inv <- solve(H)
  V_I <- U_t%*%W2%*%t(U_t)
  V_II <- U_t%*%W%*%Delta_tild%*%W%*%t(U_t)
  V_tot <- V_I+V_II
  V_II_e <- U_t%*%W%*%Delta_tild_e%*%W%*%t(U_t)
  V_II_nnc <- U_t%*%W%*%Delta_tild_c%*%W%*%t(U_t)
  ## MD and KC correction
  Var <- (H_inv%*%V_tot%*%H_inv)
  Var_df <- K_frac*(H_inv%*%V_tot%*%H_inv)
  Var_nnc <- H_inv%*%(V_I+V_II_nnc)%*%H_inv
  Var_nnc_df <- K_frac*Var_nnc
  Var_e <- H_inv%*%V_II_e%*%H_inv
  ##
  SD <- sqrt(diag(Var))
  SD_df <- sqrt(diag(Var_df))
  SD_nnc <- sqrt(diag(Var_nnc))
  SD_nnc_df <- sqrt(diag(Var_nnc_df))
  SD_MD_KC <- getMD_KC(K, n_k, D_k, D_k_t, Vk_inv, Wts_k, H_inv, V_II_e, p, Delta, Y_mu_k)
  SD_KC <- SD_MD_KC$SD_KC
  return(SD_KC)
}

getMD_KC <- function(K, n_k, D_k, D_k_t, Vk_inv, Wts_k, H_inv, V_II_e, p, Delta, Y_mu_k){
  
  F.MD.1 <- matrix(0, p, p)
  F.KC.1 <- matrix(0, p, p)
  F.MD.2 <- matrix(0, p, p)
  F.KC.2 <- matrix(0, p, p)
  
  Ident_k <- lapply(1:K, function(x) diag(n_k[[x]]))
  H_k <- lapply(1:K, function(x) 
    D_k[[x]]%*%H_inv%*%D_k_t[[x]]%*%Vk_inv[[x]]%*%(diag(n_k[[x]])*Wts_k[[x]]))
  
  M <- lapply(1:K, function(x) Ident_k[[x]] - H_k[[x]])
  Bkk_scalar <- lapply(1:K, function(x) Wts_k[[x]][1]*Delta[[x]][1,1]*Wts_k[[x]][1])
  
  for(k in 1:K){
    M.inv <- tryCatch(solve(M[[k]]), error=function(e) NA)
    ## correction to VI
    if(is.na(M.inv)[1]){
      # no correction to this component k
      U.wk1.MD <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]]) 
      F.MD.1   <- F.MD.1 + Wts_k[[k]][1]*(U.wk1.MD %*% t(U.wk1.MD))
      comp1k.MD <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]])
    } else if(abs(M.inv[1,1]) > 1e+10){
      #no correction to this component k
      U.wk1.MD <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]]) 
      F.MD.1   <- F.MD.1 + Wts_k[[k]][1]*(U.wk1.MD %*% t(U.wk1.MD))
      comp1k.MD <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]]) 
    } else{
      U.wk1.MD <-  D_k_t[[k]] %*% Vk_inv[[k]] %*% M.inv %*% matrix(Y_mu_k[[k]]) 
      F.MD.1   <- F.MD.1 + Wts_k[[k]][1]*(U.wk1.MD %*% t(U.wk1.MD))
      comp1k.MD <- D_k_t[[k]] %*% Vk_inv[[k]] %*% M.inv %*% matrix(Y_mu_k[[k]])
    }
    
    if(is.na(sqrtm(M[[k]])[1])){
      # no correction to this component k
      U.wk1.KC <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]]) 
      F.KC.1   <- F.KC.1 + Wts_k[[k]][1]*(U.wk1.KC %*% t(U.wk1.KC))
      comp1k.KC <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]]) 
    } else if(abs(sqrtm(M[[k]])[1,1]) > 1e+10){
      #no correction to this component k
      U.wk1.KC <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]])
      F.KC.1   <- F.KC.1 + Wts_k[[k]][1]*(U.wk1.KC %*% t(U.wk1.KC))
      comp1k.KC <- D_k_t[[k]] %*% Vk_inv[[k]] %*% matrix(Y_mu_k[[k]])
    } else{
      U.wk1.KC <-  D_k_t[[k]] %*% Vk_inv[[k]] %*% sqrtm(M.inv) %*% matrix(Y_mu_k[[k]])
      F.KC.1   <- F.KC.1 + Wts_k[[k]][1]*(U.wk1.KC %*% t(U.wk1.KC))
      comp1k.KC <- D_k_t[[k]] %*% Vk_inv[[k]] %*% sqrtm(M.inv) %*% matrix(Y_mu_k[[k]]) 
    }
    ##correction to VII (in two parts)
    F.MD.2 <- F.MD.2 + Bkk_scalar[[k]]*comp1k.MD%*%t(comp1k.MD)
    F.KC.2 <- F.KC.2 + Bkk_scalar[[k]]*comp1k.KC%*%t(comp1k.KC)
  }
  
  V_tot_MD <- F.MD.1 + F.MD.2 + V_II_e
  V_nnc_MD <- F.MD.1 + F.MD.2
  Var_MD <- (H_inv %*% V_tot_MD %*% H_inv)
  Var_nnc_MD <- (H_inv %*% V_nnc_MD %*% H_inv)
  SD_MD <- sqrt(diag(Var_MD))
  SD_nnc_MD <- sqrt(diag(Var_nnc_MD))
  
  V_tot_KC <- F.KC.1 + F.KC.2 + V_II_e
  V_nnc_KC <- F.KC.1 + F.KC.2
  Var_KC <- (H_inv %*% V_tot_KC %*% H_inv)
  Var_nnc_KC <- (H_inv %*% V_nnc_KC %*% H_inv)
  SD_KC <- sqrt(diag(Var_KC))
  SD_nnc_KC <- sqrt(diag(Var_nnc_KC))
  
  return(list(SD_MD=SD_MD, SD_nnc_MD=SD_nnc_MD, SD_KC=SD_KC, SD_nnc_KC=SD_nnc_KC))
}

expit <- function(x) exp(x)/(1+exp(x))
