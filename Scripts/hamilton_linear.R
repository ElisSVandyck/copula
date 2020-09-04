#Linear hamilton filter switching model

# 1) Obtain some data X[n x 1]

#state probabilities
#xi_j_t = Pr(s_t = j | information @ t, params)
#xi_j_t-1 = Pr(s_t-1 = i | information @ t - 1, params)

#densitites under regimes
#eta_j_t = f(y_t | s_t = j, information @ t-1, params)
# 		 = 1/(sqrt(2pi)sigma)*exp[-(y_t - beta0 - beta_1_j y_t-1)^2 / (2 sigma^2)]

#conditional density
#f(y_t|omega_t-1; theta) = sum^2_i=1 sum^2_j=1 p_i_j xi_i_t-1 eta_j_t
#desired update
#xi_j_t = sum^2_i=1 p_i_j xi_i_t-1 eta_j_t / f(y_t | information @ t - 1, params)
#conditional log-likelihood
#log f(y_1, y_2, ..., y_T | y_0; theta) = sum^T_t=1 log f(y_t | information @ t - 1, params)

#params_format <- c(beta0, beta1_1, sigma_1, beta1_2, sigma_2, p11_input, p22_input)

start_params <- c(rep(1,5), rep(0,2))
n <- length(y)
reg <- 2
xi <- matrix(0, n, reg)
xi[1,] <- pnorm(0)
fyt_store <- rep(0, n)

hamilton_simple_optimiser <- function(params){
	p11 <- pnorm(params[6])
	p22 <- pnorm(params[7])
	for (i in 2:n){
		eta_1 <- (1/(sqrt(2*pi)*params[3]))*exp(-((y[i]-params[0]-params[1]*y[i-1])^2)/(2*params[3]))
    	eta_2 <- (1/(sqrt(2*pi)*params[5]))*exp(-((y[i]-params[0]-params[1]*y[i-1])^2)/(2*params[5]))
    	fyt <- (p11*xi[(i-1),1]*eta_1)+((1 - p22)*xi[(i-1),2]*eta_1) + 
				((1 - p11)*xi[(i-1),1]*eta_2)+(p22*xi[(i-1),2]*eta_2)
	    xi[i,1] <- (p11*xi[(i-1),1]*eta_1+(1-p22)*xi[(i-1),2]*eta_1)/fyt
	    xi[i,2] <- ((1-p11)*xi[(i-1),1]*eta_2+(p22)*xi[(i-1),2]*eta_2)/fyt
	    fyt_store[i] <- fyt
	}
	return(-sum(log(fyt_store)))
}

estimates <- stats::optim(par = start_params, fn =hamilton_simple_optimiser, method = 'Nelder-Mead')




