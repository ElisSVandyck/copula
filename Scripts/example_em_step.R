library(tidyverse)
options(warn=-1)
library(copula)
library(VineCopula)
library(HMMcopula)
options(warn=0)
bac_raw <- read.csv('bac.csv')
c_raw <- read.csv('c.csv')
gs_raw <- read.csv('gs.csv')
ms_raw <- read.csv('ms.csv')
jpm_raw <- read.csv('jpm.csv')
returns <- rbind(bac_raw %>% select(Date, Close) %>% mutate(Name = 'bac'),
      c_raw %>% select(Date, Close) %>% mutate(Name = 'c'),
      jpm_raw %>% select(Date, Close) %>% mutate(Name = 'jpm')) %>%
    tibble() %>% 
    pivot_wider(id_cols = Date, names_from = Name, values_from = Close) %>%
    arrange(Date) %>% 
    na.omit() %>%
    mutate(bac_returns = (bac-lag(bac))/lag(bac)) %>%
    mutate(c_returns = (c-lag(c))/lag(c)) %>%
    mutate(jpm_returns = (jpm-lag(jpm))/lag(jpm)) %>%
    na.omit() 
X <- returns %>% select(x1 = bac_returns,
                        x2 = c_returns,
                        x3 = jpm_returns)
U <- X %>% filter(x1 != 0 & x2 != 0 & x3 != 0) %>% #Filtration optional ~ 10% of sample omitted
    pobs() %>% 
    data.frame() %>% tibble() %>% 
    select(u1 = x1, u2 = x2, u3 = x3)
n <- dim(U)[1]
d <- dim(U)[2]
reg <- 2
r <- reg * d
n0 <- floor(n/reg)
ind0 <- 1:n0
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
vine <- RVineCopSelect(
    data = U,
    familyset = 2,
    Matrix = structure_matrix
)
theta_start_reg1 <- c(vine$par,vine$par2) #par - t tau, par2 - t dof
theta_start_reg2 <- c(vine$par,vine$par2)

#Another school of though is to simply initialise the algorithm with
alpha_par1 <- rep(0, reg * 3)
alpha_par2 <- rep(log(5), reg * 3)
tau <- rep(0, reg * 3)
Q0 <- matrix(1, reg, reg)/reg
Q <- Q0

#EM Step Function 1 iteration example
#functions y, theta, Q, family

n <- dim(U)[1] #length of series
r <- dim(Q)[1] #number of regimes
eta_bar <- matrix(0, n, r)
eta <- matrix(0, n, r)
lambda <- matrix(0, n, r)
c = matrix(0, n, r)
Lambda <- array(0, c(r,r,n))
M <- matrix(0, r, r)

# Estimate a candidate R-Vine
# Structure vector 1
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
family_vector <- c(2,2,2,0,2,2,0,0,2)
family_matrix <- matrix(family_vector, 3, 3, byrow = FALSE)

theta_start_reg1 <- list(par_1_1 = vine$par, par_1_2 = vine$par2)
theta_start_reg2 <- list(par_2_1 = vine$par, par_2_2 = vine$par2)

vine_structure_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
           par = theta_start_reg1$par_1_1, par2 = theta_start_reg1$par_1_2)
vine_structure_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
           par = theta_start_reg2$par_2_1, par2 = theta_start_reg2$par_2_2)

c[,1] <- RVinePDF(newdata = U, RVM = vine_structure_1, verbose = TRUE)
c[,2] <- RVinePDF(newdata = U, RVM = vine_structure_2, verbose = TRUE)

eta_bar[n,] <- 1/r

v <-  (eta0 %*% Q) * c[1,] #numerator of eta at t = 1
eta[1,] <-  v/sum(v)

for (i in 2:n){
    v <-  ( eta[i-1, ] %*% Q) * c[i,] #numerator of eta at t > 1
    eta[i,] <-  v/sum(v)
}

v <- eta*eta_bar #Numerator of lambda (blending the foward backward)

sv0 <- rowSums(v)

for (j in 1:r){
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


for (j in 1:r){
    sv = rowSums(Lambda[j,,], dims=1)
    ssv = sum(sv)
    Qnew[j,] = sv/ssv
}

#Now for M -step of algorithm 

# Parameter set for regime 1
theta_1 <- c(vine_structure_1$par[2,1], #2,3 : 1 join
            vine_structure_1$par[3,1], #1,3 join
            vine_structure_1$par[3,2], #1,2 join
             vine_structure_1$par2[2,1], #dof on 2,3 :1 join
             vine_structure_1$par2[3,1], #dof on 2,3 :1 join
             vine_structure_1$par2[3,1]) #dof on 2,3 :1 join
#Parameter set for regime 2
theta_2 <- c(vine_structure_2$par[2,1], #2,3 : 1 join
            vine_structure_2$par[3,1], #1,3 join
            vine_structure_2$par[3,2], #1,2 join
             vine_structure_2$par2[2,1], #dof on 2,3 :1 join
             vine_structure_2$par2[3,1], #dof on 1,3 join
             vine_structure_2$par2[3,1]) #dof on 1,2 join

theta_comb <- c(theta_1, theta_2)




fun_opt_lite <- function(thetaa){
    par_1 <-  matrix(c(0, thetaa[1], thetaa[2], 0, 0, thetaa[3], 0, 0, 0), 3, 3)
    par_1_2 <- matrix(c(0, thetaa[4], thetaa[5], 0, 0, thetaa[6], 0, 0, 0), 3, 3)
    log_likelihood_1 = -sum(lambda[,1] * RVineLogLik(data = U, 
         RVM = vine_structure_1, 
         par = par_1,
         par2 = par_1_2)$V$direct)
    par_2 <- matrix(c(0, thetaa[7], thetaa[8], 0, 0, thetaa[9], 0, 0, 0), 3, 3)
    par_2_2 <- matrix(c(0, thetaa[10], thetaa[11], 0, 0, thetaa[12], 0, 0, 0), 3, 3)
    log_likelihood_2 = -sum(lambda[,2] * RVineLogLik(data = U, 
         RVM = vine_structure_2, 
         par = par_2,
         par2 = par_2_2)$V$direct)
    log_likelihood <- log_likelihood_1 + log_likelihood_2
    return(log_likelihood)
}


fun_opt_lite(thetaa = c(0,0,0,5,5,5,0.1,0.1,0.1,5,5,5)) #this returns correct log likelihood

# However stats::optim