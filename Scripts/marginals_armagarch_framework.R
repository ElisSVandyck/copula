### Estimation using ARMA-GARCH
library(copula)
library(rugarch)

data <- X %>% filter(x1 != 0 & x2 != 0 & x3 != 0)
plot.ts(X, plot.type = 'multiple', main = 'Returns')
#Can use xts for better plots
library(GGally)
ggpairs(data %>% select(bac = x1, c = x2, jpm = x3))

#Fit ARMA-GARCH to marginals
meanModel <- list(armaOrder = c(1,1))
varModel <- list(model = "sGARCH", garchOrder = c(1,1))

uspec <- ugarchspec(varModel, mean.model = meanModel, distribution.model = 'std') 

fit <- apply(data, 2, function(x) ugarchfit(uspec, data = x))

x1sim <- ugarchsim(fit = fit$x1, n.sim = 100)
x1sim <- cbind(x1sim@simulation$sigmaSim,
      x1sim@simulation$seriesSim,
      x1sim@simulation$residSim)
plot.ts(x1sim)
#Fit t-univaraite to residuals
#Take residuals to [0, 1] via t-cdf

#Obtain residuals via
fit$x1@fit$residuals
#Fit residuals to distribution via univariateML methodology
#See no-arma-no-marginal approach
#Use said estimated distribution to convert to [0,1]

#### COPULAS HERE #####

#For estimated copula params simulate out
#Use this simulated as custom distribution in ugarchsim
#simulate heavily as a bootstrap method
#obtain 99% value for VaR metric













