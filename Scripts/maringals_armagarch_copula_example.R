#script 5
X <- tibble(
  Date = returns$Date,
  jpm = returns$jpm,
  wfc = returns$wfc
)
Xzoo <- zoo(X[,c('jpm','wfc')], X$Date)
plot.zoo(Xzoo)
Y <- X[,c('jpm','wfc')]
Y <- na.omit(Y)
Y <- Y[apply(Y, 1, function(row) all(row != 0)),]
#7425 entries
#specify arma-garch
meanModel <- list(armaOrder = c(1,1))
varModel <- list(model = "sGARCH", garchOrder = c(1,1))
uspec <- ugarchspec(varModel, mean.model = meanModel,
                    distribution.model = 'sstd')

fits <- apply(Y, 2, function(Y) ugarchfit(uspec, data = Y))
Z <- sapply(fits, residuals, standardize = TRUE)
U <- pobs(Z)

plot(U[,1],U[,2])
pairs2(pobs(Y))
pairs2(U, labels = c('ARMAGARCH - JPM', 'WFC'))
################# Fit GUMBEL Copula for example
fit.gc <- fitCopula(gumbelCopula(dim = 2),
                    data = U, method = 'mpl')
fit.gc@estimate
gc <- fit.gc@copula
p2P(tau(gc), d = 2)
p2P(lambda(gc)['upper'], d = 2)
##########################


#script 6
#simulating paths
#Also note
nu.mar <- as.numeric(c(coef(fits$jpm)['shape'], coef(fits$wfc)['shape']))
# Use U from Y -> U and fit to T copula
fit.tc <- fitCopula(tCopula(dim = 2, dispstr = 'un'), 
                            data = U, method = 'itau.mpl')
nu <- tail(fit.tc@estimate, n = 1)
P <- p2P(head(fit.tc@estimate, n = -1))
tc <- fit.tc@copula
##########################
B <- 200
m <- ceiling(nrow(Y)/10)
Y.lst <- lapply(1:B, function(b) {
  #Simulate from the fitted copula
  U. <- rCopula(m, copula = tc)
  #Quantile-Quantile transform to standardised t distributions
  Z. <- sapply(1:2, function(j) sqrt((nu.mar[j]-2)/nu.mar[j]) * qt(U.[,j], df = nu.mar[j]))
  #Use these multivariate dependent t innovations to sample from the time series
  sim <- lapply(1:2, function(j)
    ugarchsim(fits[[j]], n.sim = m, m.sim = 1,
              custom.dist = list(name = "sample",
                                 distfit = Z.[,j,drop=FALSE]))
    )
  sapply(sim, function(x) fitted(x))
})
Y.matrix <- array(as.numeric(unlist(Y.lst)), dim = c(m,2,200))

library(RColorBrewer)

plot.ts(Y.matrix[,1,1:10], plot.type = 'single', col = brewer.pal(n = 10, name = "Greens"))


#script 7
Ys <- rowSums(Y)
Ys. <- sapply(Y.lst, rowSums)
Ys.mean <- rowMeans(Ys.)
Ys.CI <- apply(Ys., 1, function(x) quantile(x, probs = c(0.025, 0.975)))
alpha <- 0.99
VaR <- apply(Ys., 1, function(x) quantile(x, probs = alpha))
#plot
n <- nrow(Y)
tm <- 1:(start+n+m-1)
start <- tm[1]
past <- tm[start:(start+n-1)]
future <- tm[(start+n):(start+n+m-1)]


plot(past, Ys, type = "l", xlim = range(c(past, future)), xlab = "", ylab = "") # actual (past) losses
polygon(c(future, rev(future)), c(Ys.CI[1,], rev(Ys.CI[2,])),
        border = NA, col = "grey80") # CI region
lines(future, Ys.mean, col = "royalblue3") # predicted aggregated loss
lines(future, Ys.CI[1,], col = "grey50") # lower CI
lines(future, Ys.CI[2,], col = "grey50") # upper CI
lines(future, VaR, col = "maroon3") # VaR_alpha
legend("bottomright", bty = "n", lty = rep(1, 4),
       col = c("black", "royalblue3", "grey50", "maroon3"),
       legend = c("(Aggregated) loss", "(Simulated) predicted loss",
                  "95% CIs", as.expression(substitute("Simulated"~VaR[a], list(a = alpha)))))



