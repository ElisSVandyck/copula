# Estimate n-regime bivariate t-copula by EM
# Feed in U generated by import script
# Some print statements are added to the code to speed up trouble shooting,
# feel free to delete.

# Currently only works with reg = 2

tcop_em <- function(U, reg, max_iter){
  eps <- 0.001
  iter = 10
  n = dim(U)[1]
  d = dim(U)[2]
  U = floor(apply(U,2,rank))/(n+1)
  r = reg*d
  n0 = floor((n/reg))
  ind0 = 1:n0
  alpha0 = rep(0, reg)
  tau = rep(0, reg)
  for (j in 1:reg){
    ind = (j-1)*n0 + ind0
    x = U[ind,]
    tau[j] = stats::cor(x, method = 'kendall')[1,2]
    if (tau[j] <= 0.1){
      tau[j] = 0.1
    } else if (tau[j] >= 0.9){
      tau[j] = 0.9
    }
  }
  alpha0 = ParamTau('t', tau = tau)
  alpha0 = c(alpha0, log(5))
  Q0 = matrix(1, reg, reg)/reg
  
  #Warm up params
  for (k in 1:iter){
    # Apply EM algorithm
    emstep=em_algo_step(U=U,theta=alpha0,Q=Q0)
    nu=emstep$nu
    alpha_new = emstep$theta_new
    Qnew = emstep$Qnew
    eta = emstep$eta
    eta_bar =emstep$eta_bar
    lambda= emstep$lambda
    Lambda=emstep$Lambda
    Q0 = Qnew
    alpha0 = alpha_new
    print(Qnew)
  }
  #Iterate untill likelihood gain decreased sufficiently or maxiter hit
  for (k in 1:max_iter){
    emstep= em_algo_step(U=U,theta=alpha0,Q=Q0)
    nu = emstep$nu
    alpha_new = emstep$theta_new
    Qnew = emstep$Qnew
    eta = emstep$eta
    eta_bar = emstep$eta_bar
    lambda= emstep$lambda
    Lambda = emstep$Lambda
    sum1 = sum(abs(alpha0))
    sum2 = sum(abs(alpha_new-alpha0))
    if (sum2 < sum1 * r * eps){
      break
    }
    Q0 = Qnew
    alpha0 = alpha_new
    print(Qnew)
  }
  alpha = alpha_new
  Q = Qnew
  tau = KendallTau('t', alpha)
  theta = ParamCop('t', alpha)
  # Output
  out = list(alpha = alpha, Q = Q, tau = tau, theta = theta)
  return(out)
  
}

em_algo_step <- function(U, theta, Q){
  print('Begininning EM iteration')
  print('E-Step')
  n = dim(U)[1]
  r = dim(Q)[2]
  eta_bar = matrix(0, n, r)
  eta = matrix(0, n, r)
  lambda = matrix(0, n, r)
  c = matrix(0, n, r)
  Lambda = array(0, c(r, r, n))
  M = matrix(0, r, r)
  #### Populate
  for (j in 1:r){
    c[,j] = copulaFamiliesPDF('t', U, theta[j], theta[r+1])
  }
  # Expectation Step
  # Backward
  eta_bar[n,] = 1/r
  for (k in 1:(n-1)){
    i = n-k
    j = i+1
    v = (eta_bar[j,] * c[j,]) %*% t(Q)
    eta_bar[i,] = v/sum(v)
  }
  eta0 = rep(1, r)/r
  v = (eta0 %*% Q) * c[1,]
  #Forward
  eta[1,] = v/sum(v)
  for (i in 2:n){
    v = (eta[(i-1),] %*% Q) * c[i,]
    eta[i,] = v/sum(v)
  }
  v = eta*eta_bar
  sv0 = rowSums(v)
  for(j in 1:r){
    lambda[,j] = v[,j]/sv0
  }
  gc = eta_bar * c
  M = Q * (as.matrix(eta0)%*% gc[1,] )
  MM = sum(M)
  Lambda[,,1] = M/MM
  for(i in 2:n){
    M = Q * (as.matrix(eta[i-1,]) %*% gc[i,])
    MM = sum(M)
    Lambda[,,i] = M/MM
  }
  nu = colMeans(lambda)
  Qnew = Q
  for (j in 1:r){
    sv = rowSums(Lambda[j,,], dims = 1)
    ssv = sum(sv)
    Qnew[j,] = sv/ssv
  }
  #Maximisation Step
  print('M-Step')
 # print(sum(lambda[,1]*log(copulaFamiliesPDF('t', U, theta[1], theta[r+1]))))
  print('Beginning Optimisation')
  maximisation_for_em_t_bi <- function(theta){
    log_likelihood_1 = -sum(lambda[,1]*log(copulaFamiliesPDF('t', U ,theta[1], theta[r+1])))
    log_likelihood_2 = (-sum(lambda[,2]* log(copulaFamiliesPDF('t', U, theta[2], theta[r+1]) )))
    log_likelihood = log_likelihood_1 + log_likelihood_2
    return(log_likelihood)
  }
  theta_new = stats::optim(par = theta, fn = maximisation_for_em_t_bi, method = 'Nelder-Mead')$par
  
  out = list(nu=nu, theta_new=theta_new, Qnew=Qnew, eta=eta,eta_bar=eta_bar,
             lambda=lambda, Lambda=Lambda)
  
  return(out)
}



estimates <- tcop_em(U = U, reg = 2, max_iter = 100)
