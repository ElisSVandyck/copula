U <- X %>% filter(x1 != 0 & x2 != 0 & x3 != 0) %>% #Filtration optional ~ 10% of sample omitted
  pobs() %>% 
  data.frame() %>% tibble() %>% 
  select(u1 = x1, u2 = x2, u3 = x3)


n <- dim(U)[1]
d <- dim(U)[2]
reg <- 2
r <- reg * d
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
vine <- RVineCopSelect(
  data = U,
  familyset = 2,
  Matrix = structure_matrix
)
theta_1_pre <- c(vine$par[2,1],vine$par[3,1],vine$par[3,2]) #par - t tau, par2 - t dof
theta_2_pre <- c(theta_1_pre)*0.8
dof <- log(10)
theta_comb <- c(theta_1_pre, theta_2_pre)
Q0 = matrix(1,reg,reg)/reg
theta_comb0 <- theta_comb
# Classic EM Step inputs are:
# y, theta, Q
# EMStep(U, theta_comb, Q0)
U <- U; theta_comb <- theta_comb0; Q <- Q0
#EM Step function
input_par_1 <- matrix(c(0, theta_comb[1], theta_comb[2], 0, 0, 
                        theta_comb[3], 0, 0, 0), 3, 3)
input_par_2 <- matrix(c(0, theta_comb[4], theta_comb[5], 0, 0, 
                        theta_comb[6], 0, 0, 0), 3, 3)
input_dof <- matrix(c(0, dof, dof, 0, 0, 
                      dof, 0, 0, 0), 3, 3)
vine_structure_update_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                       par = input_par_1, par2 = input_dof)
vine_structure_update_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                       par = input_par_2, par2 = input_dof)
#Input checks completed
n = dim(U)[1]
r = dim(Q)[2]
eta_bar = matrix(0, n, r)
eta = matrix(0, n, r)
lambda = matrix(0, n, r)
c = matrix(0, n, r)
Lambda = array(0, c(r, r, n))
M = matrix(0, r, r)
###########
c[,1] = RVinePDF(newdata = U, RVM = vine_structure_update_1, verbose = TRUE)
c[,2] = RVinePDF(newdata = U, RVM = vine_structure_update_2, verbose = TRUE)
#E Step
eta_bar[n,] <- 1/r
for (k in 1:(n-1)){
  i = n-k
  j = i+1
  v = (eta_bar[j,] * c[j,]) %*% t(Q)
  eta_bar[i,] = v/sum(v)
}
eta0 = rep(1,r)/r
v <-  (eta0 %*% Q) * c[1,] #numerator of eta at t = 1
eta[1,] <-  v/sum(v)
for (i in 2:n){
  v <-  ( eta[i-1, ] %*% Q) * c[i,] #numerator of eta at t > 1
  eta[i,] <-  v/sum(v)
}

v <- eta*eta_bar #Numerator of lambda (blending the foward backward)
sv0 <- rowSums(v)
for (j in 1:2){
  lambda[,j] <- v[,j]/sv0
}

gc <- eta_bar * c

M <- Q * (as.matrix(eta0) %*% gc[1,])

MM <- sum(M)

Lambda[,,1] <- M/MM

for (i in 2:n){
  M = Q * ( as.matrix(eta[i-1,]) %*% gc[i,]) #Numerator of lambda
  MM = sum(M)
  Lambda[,,i] = M/MM 
}

nu <- colMeans(lambda)
Qnew <- Q 


for (j in 1:2){
  sv = rowSums(Lambda[j,,], dims=1)
  ssv = sum(sv)
  Qnew[j,] = sv/ssv
}

rvine_t_optim <- function(theta_comb){
  input_par_1 <- matrix(c(0, theta_comb[1], theta_comb[2], 0, 0, 
                          theta_comb[3], 0, 0, 0), 3, 3)
  input_par_2 <- matrix(c(0, theta_comb[4], theta_comb[5], 0, 0, 
                          theta_comb[6], 0, 0, 0), 3, 3)
  input_dof <- matrix(c(0, dof, dof, 0, 0, 
                        dof, 0, 0, 0), 3, 3)
  vine_structure_update_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                         par = input_par_1, par2 = input_dof)
  vine_structure_update_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                         par = input_par_2, par2 = input_dof)
  loglik_1 <- RVinePDF(newdata = U, RVM = vine_structure_update_1)
  loglik_1sum <- -sum(lambda[,1]*log(loglik_1))
  loglik_2 <- RVinePDF(newdata = U, RVM = vine_structure_update_2)
  loglik_2sum <- -sum(lambda[,2]*log(loglik_2))
  return((loglik_1sum + loglik_2sum))
}
results <- stats::optim(par = theta_comb, fn = rvine_t_optim, 
                        method = 'Nelder-Mead',
                        control = list(maxit = 1000))

