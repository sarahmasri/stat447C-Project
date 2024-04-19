library(rstan)
library(glmnet)
library(MASS)

set.seed(447)
setwd("/Users/sarahmasri/Desktop/STAT  447C/stat447C-Project/Code/")

# Generate synthetic data
n <- 1000  # Number of samples
p <- 100   # Number of predictors

# True coefficients
beta_true <- numeric(p)
ind <- sample(1:100, 7)
beta_true[ind] <- runif(7, min=-5, max=5)

# Generate predictors
X <- matrix(rnorm(n * p), nrow = n)

# Generate response variable
y <- X %*% beta_true + rnorm(n, sd = 0.5)

# Normalize predictors
X <- scale(X)

# Elastic Net
lambda_seq <- 10^seq(3, -3, length.out = 100)  # Sequence of lambda values
EN_fit <- cv.glmnet(X, y, alpha = 0.5, lambda = lambda_seq)

# Extract optimal lambda
lambda_opt_EN <- EN_fit$lambda.min

# Plot most optimal lambda for elastic net
plot(EN_fit)
abline(v = log(lambda_opt_EN), col = "red", lty = 2)
legend("bottomright", legend = "Optimal Lambda", col = "red", lty = 2, yjust=-1) 


# Fit final elastic net model with optimal lambda
EN_model <- glmnet(X, y, alpha = 0.5, lambda = lambda_opt_EN)





# Compile Stan model
file <- "./horseshoe.stan"  # Paste the Stan model code here
stan_model <- stan_model(file = file)


# Prepare data for Stan
stan_data <- list(
  n = n,
  p = p,
  X = X,
  y = as.vector(y)
)


# Fit Bayesian regression model
bayes_fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4)


# Predictive performance
EN_pred <- predict(EN_model, newx = X)
bayes_pred <- extract(bayes_fit)$y_pred

# Mean squared error
EN_mse <- mean((y - EN_pred)^2)
bayes_mse <- mean((y - colMeans(bayes_pred))^2)

# Variable selection
EN_selected <- coef(EN_model) != 0

thresh <- 0.001
bayes_coefs <- as.vector(summary(bayes_fit, pars = c("beta"))$summary[,1])
bayes_selected <- bayes_coefs > thresh

par(mfrow=c(1, 2))
hist(as.vector(coef(EN_model)), col = 'skyblue3', breaks=20, xlab = 'beta', main="NEN Coefficients")
hist(extract(bayes_fit)$beta, col = 'skyblue3', breaks=20,  xlab = 'beta', main="Horseshoe Coefficients")



# Print results
cat("Elastic Net MSE:", EN_mse, "\n")
cat("Bayesian Regression MSE:", mean(bayes_mse), "\n\n")
cat("Elastic Net Selected Variables:", sum(EN_selected), "\n")
cat("Bayesian Regression Selected Variables:", sum(bayes_selected), "\n")
cat("Elastic Net Correctly Selected Variables:", sum(which(EN_selected) %in% ind ), "\n")
cat("Bayesian Regression Correctly Selected Variables:", sum(which(bayes_selected) %in% ind ), "\n")
