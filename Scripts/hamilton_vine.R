library(tidyverse)
library(copula)
library(VineCopula)
library(HMMcopula)
bac_raw <- read.csv('bac.csv')
c_raw <- read.csv('c.csv')
gs_raw <- read.csv('gs.csv')
ms_raw <- read.csv('ms.csv')
jpm_raw <- read.csv('jpm.csv')
returns <- rbind(bac_raw %>% select(Date, Close) %>% mutate(Name = 'bac'),
                 c_raw %>% select(Date, Close) %>% mutate(Name = 'c'),
                 jpm_raw %>% select(Date, Close) %>% mutate(Name = 'jpm'))%>%
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


uu <- X %>% pobs() %>% 
  data.frame() %>% tibble() %>% 
  select(u1 = x1, u2 = x2, u3 = x3)
u <- as.matrix(uu[6000:dim(uu)[1],])  # here we are just reducing the number of obs

#theta intialise
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
vine <- RVineCopSelect(
  data = u,
  familyset = 2,
  Matrix = structure_matrix
)

theta_start_reg1 <-  list(par_1_1 = (vine$par), par_1_2 = (vine$par2))
theta_start_reg2 <- list(par_2_1 = (vine$par*0.4), par_2_2 = (vine$par2*1.2))
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
family_vector <- c(2,2,2,0,2,2,0,0,2)
family_matrix <- matrix(family_vector, 3, 3, byrow = FALSE)
vine_structure_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                par = theta_start_reg1$par_1_1, par2 = theta_start_reg1$par_1_2)
vine_structure_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                par = theta_start_reg2$par_2_1, par2 = theta_start_reg2$par_2_2)
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
p11 <- 0.5
p22 <- 0.5
theta_plus_P <- c(theta_comb, p11, p22)


input_par_1_1 <-  matrix(c(0, theta_plus_P[1], theta_plus_P[2], 0, 0, 
                           theta_plus_P[3], 0, 0, 0), 3, 3)
input_par_1_2 <- matrix(c(0, theta_plus_P[4], theta_plus_P[5], 0, 0, 
                          theta_plus_P[6], 0, 0, 0), 3, 3)
input_par_2_1 <-  matrix(c(0, theta_plus_P[7], theta_plus_P[8], 0, 0, 
                           theta_plus_P[9], 0, 0, 0), 3, 3)
input_par_2_2 <- matrix(c(0, theta_plus_P[10], theta_plus_P[11], 0, 0, 
                          theta_plus_P[12], 0, 0, 0), 3, 3)
vine_structure_update_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                       par = input_par_1_1, par2 = input_par_1_2)
vine_structure_update_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                       par = input_par_2_1, par2 = input_par_2_2)
#starting parameter values
param0=c(theta_comb,0,0)

reg <- 2
n <- dim(u)[1]
u <- as.matrix(u)  # here we are just reducing the number of obs
xi <- matrix(0, n, reg)
xi[1,] <- pnorm(0)
fyt_store <- rep(0, n)

hamilton_filter_map <- function(param){
  #Update parameters
  par_11 <-  matrix(c(0, param[1], param[2], 0, 0, 
                             param[3], 0, 0, 0), 3, 3)
  par_12 <- matrix(c(0, param[4], param[5], 0, 0, 
                            param[6], 0, 0, 0), 3, 3)
  par_21 <-  matrix(c(0, param[7], param[8], 0, 0, 
                             param[9], 0, 0, 0), 3, 3)
  par_22 <- matrix(c(0, param[10], param[11], 0, 0, 
                            param[12], 0, 0, 0), 3, 3)
  vine_structure_update_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                         par = par_11, par2 = par_12)
  vine_structure_update_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
                                         par = par_21, par2 = par_22)
  p11 <- pnorm(param[13])
  p22 <- pnorm(param[14])
  for(i in 2:n){
    #Densities under regimes
    eta_1 <- RVinePDF(newdata = u[i,], RVM = vine_structure_update_1, verbose = TRUE)
    eta_2 <- RVinePDF(newdata = u[i,], RVM = vine_structure_update_2, verbose = TRUE)
    fyt <- (p11*xi[(i-1),1]*eta_1)+((1 - p22)*xi[(i-1),2]*eta_1) + 
      ((1 - p11)*xi[(i-1),1]*eta_2)+(p22*xi[(i-1),2]*eta_2)
    #Conditional density
    xi[i,1] <- (p11*xi[(i-1),1]*eta_1+(1-p22)*xi[(i-1),2]*eta_1)/fyt
    xi[i,2] <- ((1-p11)*xi[(i-1),1]*eta_2+(p22)*xi[(i-1),2]*eta_2)/fyt
    fyt_store[i] <- fyt
  }
  
  print(param)
  print(-sum(log(fyt_store[2:n])))
  #return(-sum(log(fyt_store[2:n])))
}

time01 <- proc.time()
optim(par = param0, fn = hamilton_filter_map, method = 'L-BFGS-B',
      lower = c(rep(-0.99,3),rep(2.5,3),rep(-0.99,3),rep(2.5,3),rep(-Inf,2)),
      upper = c(rep(0.99,3),rep(100,3),rep(0.99,3),rep(100,3),rep(Inf,2)),
      control = list(maxix = 100))
time11<-proc.time() - time01



