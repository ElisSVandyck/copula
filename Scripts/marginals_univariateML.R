# Fitting marginals without a measure of univariate dependence (lag based dependence)

install.packages('univariateML')
library(univariateML)
data <- X %>% 
  filter(x1 != 0 & x2 !=0 & x3 != 0) %>% 
  data.frame()
nrow(data)
margin_fits <- lapply(data, model_select, models = univariateML_models, criterion = 'aic')
U. <- lapply(seq_along(data), function(j) pml(data[[j]], margin_fits[[j]]))
U <- data.frame(u1 = U.[[1]], u2 = U.[[2]], u3 = U.[[3]])
#theta intialise
structure_vector <- c(3,2,1,0,2,1,0,0,1)
structure_matrix <- matrix(structure_vector, 3, 3, byrow = FALSE)
vine <- RVineCopSelect( 
  data = U,
  familyset = 2,
  Matrix = structure_matrix
)
#RVine selection criteria needs to be a chapter
#plot(vine)

Usim <- as.data.frame(RVineSim(N = nrow(U), RVM = vine))
Usim_out <- lapply(seq_along(Usim), function(j) qml(data.frame(Usim)[[j]], margin_fits[[j]]))
Usim_out <- data.frame(u1sim = Usim_out[[1]],u2sim = Usim_out[[2]],u3sim = Usim_out[[3]])

pairs(Usim_out)
