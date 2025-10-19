

#### Covariance matrix ####
stbl_param_covmtrx <- function(t_par, theta){

  phi_2t <- ComplexCF(t=2*t_par, theta=theta, pm=1)
  phi_m2t <- ComplexCF(t=-2*t_par, theta=theta, pm=1)

  phi_t <- ComplexCF(t=t_par, theta=theta, pm=1)
  phi_mt <- ComplexCF(t=-1*t_par, theta=theta, pm=1)

  el_11 <- Re(1/4 * (phi_2t + 2 + phi_m2t - (phi_t)^2 - 2*phi_t*phi_mt - phi_mt^2))
  el_22 <- Re(1/4 * ((phi_t)^2 - 2*phi_t*phi_mt + (phi_mt)^2) - 1/4*(phi_2t + phi_m2t - 2))
  el_nd <- Re(1/(4i) * (phi_2t - (phi_t)^2 - phi_m2t + (phi_mt)^2))

  matrix(c(el_11, el_nd, el_nd, el_22), 2, 2)
}


# t_par=1; theta <- c(1.5, 1, 3, 10)
# mtrx <- stbl_param_covmtrx(t_par=4, theta)
# solve(mtrx)


Z_0 <- function(t_par=4, theta){
  phi_t <- ComplexCF_cpp(t=t_par, theta=theta)
  c(Re(phi_t), Im(phi_t))
}


Z_n <- function(t_par, X){
  phi_t <- sum(exp(1i*t_par*X))/length(X)
  c(Re(phi_t), Im(phi_t))
}


############# Abandoned for the moment ####################
# Not used ??
ECF_Stab_Distance <- function(X, t_par, theta){

  Z_n <- Z_n(t_par=t_par, X=X)
  Z_0 <- Z_0(t_par=t_par, theta)
  Sigma <- stbl_param_covmtrx(t_par=t_par, theta)
  as.numeric(t(Z_n - Z_0) %*% solve(Sigma) %*% (Z_n - Z_0))
}

# X <- rstable(n=1e5, alpha=theta[1], beta=theta[2], gamma=theta[3], delta=theta[4])
# ECF_Stab_Distance(X, t_par, theta)

# X <- rstable(n=1e5, alpha=1.1, beta=0*theta[2], gamma=theta[3], delta=theta[4])
# ECF_Stab_Distance(X, t_par, theta)


####
# Not used ??
#
# library(deforeStable)
# data("geoimages")
#
# obj <- geoimages[[2]]
# plotRGB(obj, scale=1, asp=1)
#
# mtrx <- as.matrix(obj)
# pars <- Koutparams(mtrx)
#
# dd_cvm <- Forest_Tester_stab_ECF(params=pars, dataset=mtrx, t_par=1)
# dd_cvm
Forest_Tester_stab_ECF <- function(params, dataset, t_par){

  cvm_red <- ECF_Stab_Distance(X=as.vector(dataset[,1]), t_par, as.vector(params['red',]))
  cvm_green <- ECF_Stab_Distance(X=as.vector(dataset[,2]), t_par, as.vector(params['green',]))
  cvm_blue <- ECF_Stab_Distance(X=as.vector(dataset[,3]), t_par, as.vector(params['blue',]))

  data.frame(t(c(cvm_red, cvm_green, cvm_blue)), fix.empty.names = FALSE)
}


#####
# Not used ??
Multiple_ECF_Tester <- function(data, params, t_par){

  test <- plyr::ldply(params, Forest_Tester_stab_ECF, dataset = data, t_par=t_par,
                      .parallel = FALSE, .paropts = list(.packages = 'stabledist'))

  test <- cbind(test, sum = rowSums(test))
  res <- test[which.min(test$sum),]
  res
}
