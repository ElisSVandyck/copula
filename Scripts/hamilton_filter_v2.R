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
    #theta intialise
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
vine <- RVineCopSelect(
    data = U,
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
reg <- 2
n <- dim(U)[1]
xi <- matrix(0, n, reg)
p11 <- 0.6
p22 <- 0.4
U <- as.matrix(U)
theta_plus_P <- c(theta_comb, qnorm(p11), qnorm(p22))

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

hamilton_filter_map <- function(theta_plus_P_input){
    print('Beginning function')
    xi <- matrix(0, n, reg)
    xi[1,] <- 0.5
    fyt_store <- rep(0, n)
        #Update parameters
    input_par_1_1 <-  matrix(c(0, theta_plus_P_input[1], theta_plus_P_input[2], 0, 0, 
                           theta_plus_P_input[3], 0, 0, 0), 3, 3)
    input_par_1_2 <- matrix(c(0, theta_plus_P_input[4], theta_plus_P_input[5], 0, 0, 
                              theta_plus_P_input[6], 0, 0, 0), 3, 3)
    input_par_2_1 <-  matrix(c(0, theta_plus_P_input[7], theta_plus_P_input[8], 0, 0, 
                               theta_plus_P_input[9], 0, 0, 0), 3, 3)
    input_par_2_2 <- matrix(c(0, theta_plus_P_input[10], theta_plus_P_input[11], 0, 0, 
                              theta_plus_P[12], 0, 0, 0), 3, 3)
    vine_structure_update_1 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
               par = input_par_1_1, par2 = input_par_1_2)
    vine_structure_update_2 <- RVineMatrix(Matrix = structure_matrix, family = family_matrix, 
               par = input_par_2_1, par2 = input_par_2_2)
    p11 <- pnorm(theta_plus_P_input[13])
    p22 <- pnorm(theta_plus_P_input[14])
    for(i in 2:n){
    #Densities under regimes
        if (i %% 1000 == 0){
            print(i)
        }
        eta_1 <- RVinePDF(newdata = U[i,], RVM = vine_structure_update_1, verbose = TRUE)
        eta_2 <- RVinePDF(newdata = U[i,], RVM = vine_structure_update_2, verbose = TRUE)
        fyt <- (p11 * xi[(i-1),1] * eta_1)+((1 - p22) * xi[(i-1),2] * eta_1) + 
            ((1 - p11 ) * xi[(i-1),1] * eta_2)+(p22 * xi[(i-1),2] * eta_2)
        #Conditional density
        xi[i,1] <- (p11*xi[(i-1),1]*eta_1+(1-p22)*xi[(i-1),2]*eta_1)/fyt
        xi[i,2] <- ((1-p11)*xi[(i-1),1]*eta_2+(p22)*xi[(i-1),2]*eta_2)/fyt
        fyt_store[i] <- fyt
    }
    print(theta_plus_P_input)
    print(-sum(log(fyt_store[2:n])))
    #log_fyt_store <- log(fyt_store)
    return(-sum(log(fyt_store[2:n])))
}


stats::optim(par = theta_plus_P, fn = hamilton_filter_map, method = 'Nelder-Mead',
            control = list(maxit = 20, trace = TRUE))