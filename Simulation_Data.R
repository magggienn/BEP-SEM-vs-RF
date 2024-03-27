
### Monte Carlo Simulation for CCLSLV method for prediction 
### ReMa IDA Traineeship III
### Author: Anouk Bouma
### Supervisor: dr. Katrijn van Deun
### Last modified: November 5th, 2023

# Loading the required functions
#source("~/Simulation/CCLSLV.R")

library(MASS)

# General parameters for every dataset
n <- c(300,30,20) # Sample sizes in each setup
J <- 20 # Times in the replication function
R <- 2 # Number of latent variables

b1 <- .5 # True value for the slope of predictor 1
b2 <- .3 # True value for the slope of predictor 2

reps <- 1000 # Number of repititions

# Creating the factor structure
p1g1 <- c(0,rep(sqrt(.6),10),rep(0,9)) # predictor 1 in group 1, a zero, ten loadings and nine zeroes
p2g1 <- c(sqrt(.6),rep(0,10),rep(sqrt(.6),9)) # predictor 2 in group 1, a loading, ten zeroes and nine loadings
P1 <- cbind(p1g1,p2g1) # make it a dataframe

# With random noise:
PSI1 <- diag(1-rowSums(P1^2)) #factors have variance 1, common+unique; this obtains unique variance
SIGMA1 <- P1%*%t(P1)+ PSI1 # Sigma with noise

# Population data without noise:
SIGMApop <- P1%*%t(P1) #Population Sigma without noise

# Empty matrix to store the performance measures
perf.meas <- matrix(NA, nrow = 3, ncol = 5)

#------------------------------------------------------#

# Setup 1 
# n > p

#------------------------------------------------------#

# Systematic component #
set.seed(23465) # Set the seed for reproducible results
par.est1 <- matrix(NA, nrow = reps, ncol = 3) # Empty matrix to store estimates

for (i in 1:reps) { # Start the loop
  j <- 1 #choose sample size condition in vector n
  
  # Population data without noise:
  datapop1 <- mvrnorm(n = n[j], mu = rep(0,J), Sigma = SIGMApop, empirical = FALSE)
  fpop1 <- CCLSLV(datapop1, R, 20, 100, .001)
  
  Ypop1 <- b1*fpop1$scores[,1] + b2*fpop1$scores[,2] # Population scores for Y
  
  # Data with random noise:
  data1 <- mvrnorm(n = n[j], mu = rep(0,J), Sigma = SIGMA1, empirical = TRUE) #mvrnorm Random noise
  
  # Estimate factor scores here: code Tra CCLSLV
  fhat1 <- CCLSLV(data1, R, 20, 100, .001) # R = 2, number of latent factors
  
  Yhat1 <- b1*fhat1$scores[, 1] + b2*fhat1$scores[, 2] + rnorm(n[j], 0, 1) # True DGP, with N(0, 1) error
  model1 <- lm(Yhat1[1:200 ] ~ fhat1$scores[1:200, 1] + fhat1$scores[1:200, 2] - 1) 
  
  vcv1 <- vcov(model1) # Variance-covariance matrix
  
  Ypred1 <- model1$coef[1]*fhat1$scores[201:300, 1] + model1$coef[2]*fhat1$scores[201:300, 2] #Prediction in testing set
  
  par.est1[i, 1] <- model1$coef[1] # Put the estimate for the first coefficient in the first column
  par.est1[i, 2] <- model1$coef[2] # Put the estimate for the coefficient on b2 in the second column
  par.est1[i, 3] <- ((sum(Ypop1[201:300] - Ypred1)^2)/(sum(Ypop1[201:300]^2))) # Prediction error
} # End the loop

# Bias -> average value of the estimate 
# Absolute Bias
ab.beta1.1 <- abs(b1 - (mean(par.est1[ ,1]))) 
perf.meas[1,1] <- ab.beta1.1

ab.beta2.1 <- abs(b2 - (mean(par.est1[ ,2]))) # Mean of the absolute difference between the estimate and the true beta for each sample
perf.meas[1,2] <- ab.beta2.1

# MSE (efficiency)
mse.beta1.1 <- mean((par.est1[ , 1] - b1) ^2) # Mean squared error
perf.meas[1,3] <- mse.beta1.1

mse.beta2.1 <- mean((par.est1[ , 2] - b2) ^2) # Mean squared error
perf.meas[1,4] <- mse.beta2.1

# Prediction error
pred.error1 <- mean(par.est1[ , 3])
perf.meas[1,5] <- pred.error1


#------------------------------------------------------#

# Setup 2
# n = p

#------------------------------------------------------#

# Systematic component #
set.seed(23465) # Set the seed for reproducible results
par.est2 <- matrix(NA, nrow = reps, ncol = 3) # Empty matrix to store estimates

for (i in 1:reps) { # Start the loop
  j <- 2 #choose sample size condition in vector n
  
  # Population data without noise:
  datapop2 <- mvrnorm(n = n[j], mu = rep(0,J), Sigma = SIGMApop, empirical = FALSE)
  fpop2 <- CCLSLV(datapop2, R, 20, 100, .001)
  
  Ypop2 <- b1*fpop2$scores[,1] + b2*fpop2$scores[,2] # Population scores for Y
  
  # Data with random noise:
  data2 <- mvrnorm(n = n[j], mu = rep(0,J), Sigma = SIGMA1, empirical = TRUE) #mvrnorm Random noise 
  
  # Estimate factor scores here: code Tra CCLSLV
  fhat2 <- CCLSLV(data2, R, 20, 100, .001) # R = 2, given two latent variables
  
  Yhat2 <- b1*fhat2$scores[, 1] + b2*fhat2$scores[, 2] + rnorm(n[j], 0, 1) # True DGP, with N(0, 1) error
  model2 <- lm(Yhat2[1:20 ] ~ fhat2$scores[1:20, 1] + fhat2$scores[1:20, 2] - 1) 
  
  vcv2 <- vcov(model2) # Variance-covariance matrix
  
  Ypred2 <- b1*fhat2$scores[21:30, 1] + b2*fhat2$scores[21:30, 2] #Prediction in testing set
  
  par.est2[i, 1] <- model2$coef[1] # Put the estimate for the first coefficient in the first column
  par.est2[i, 2] <- model2$coef[2] # Put the estimate for the coefficient on b2 in the second column
  par.est2[i, 3] <- ((sum(Ypop2[21:30] - Ypred2)^2)/(sum(Ypop2[21:30]^2))) # Prediction error
} # End the loop

# Bias -> average value of the estimate 
# Absolute Bias
ab.beta1.2 <- abs(b1 - (mean(par.est2[ ,1]))) 
perf.meas[2,1] <- ab.beta1.2

ab.beta2.2 <- abs(b2 - (mean(par.est2[ ,2]))) # Mean of the absolute difference between the estimate and the true beta for each sample
perf.meas[2,2] <- ab.beta2.2

# MSE (efficiency)
mse.beta1.2 <- mean((par.est2[ , 1] - b1) ^2) # Mean squared error
perf.meas[2,3] <- mse.beta1.2

mse.beta2.2 <- mean((par.est2[ , 2] - b2) ^2) # Mean squared error
perf.meas[2,4] <- mse.beta2.2

# Prediction error
pred.error2 <- mean(par.est2[ , 3])
perf.meas[2,5] <- pred.error2

#------------------------------------------------------#

# Setup 3
# n < p

#------------------------------------------------------#

# Systematic component #
set.seed(23465) # Set the seed for reproducible results
par.est3 <- matrix(NA, nrow = reps, ncol = 3) # Empty matrix to store estimates

for (i in 1:reps) { # Start the loop
  j <- 3 #choose sample size condition in vector n
  
  # Population data without noise:
  datapop3 <- mvrnorm(n = n[j], mu = rep(0,J), Sigma = SIGMApop, empirical = FALSE)
  fpop3 <- CCLSLV(datapop3, R, 20, 100, .001)
  
  Ypop3 <- b1*fpop3$scores[,1] + b2*fpop3$scores[,2] # Population scores of Y
  
  # Data with random noise:
  data3 <- mvrnorm(n = n[j], mu = rep(0,J), Sigma = SIGMA1, empirical = TRUE) #mvrnorm Random noise
  
  # Estimate factor scores here: code Tra CCLSLV
  fhat3 <- CCLSLV(data3, R, 20, 100, .001) # R = 2, given two latent variables
  
  Yhat3 <- b1*fhat3$scores[, 1] + b2*fhat3$scores[, 2] + rnorm(n[j], 0, 1) # True DGP, with N(0, 1) error
  model3 <- lm(Yhat3[1:15 ] ~ fhat3$scores[1:15, 1] + fhat3$scores[1:15, 2] - 1) 
  
  vcv3 <- vcov(model3) # Variance-covariance matrix
  
  Ypred3 <- b1*fhat3$scores[16:20, 1] + b2*fhat3$scores[16:20, 2] #Prediction in testing set
  
  par.est3[i, 1] <- model3$coef[1] # Put the estimate for the first coefficient in the first column
  par.est3[i, 2] <- model3$coef[2] # Put the estimate for the coefficient on b2 in the second column
  par.est3[i, 3] <- ((sum(Ypop3[16:20] - Ypred3)^2)/(sum(Ypop3[16:20]^2)))
} # End the loop

# Bias -> average value of the estimate 
# Absolute Bias
ab.beta1.3 <- abs(b1 - (mean(par.est3[ ,1]))) 
perf.meas[3,1] <- ab.beta1.3

ab.beta2.3 <- abs(b2 - (mean(par.est3[ ,2]))) # Mean of the absolute difference between the estimate and the true beta for each sample
perf.meas[3,2] <- ab.beta2.3

# MSE (efficiency)
mse.beta1.3 <- mean((par.est3[ , 1] - b1) ^2) # Mean squared error
perf.meas[3,3] <- mse.beta1.3

mse.beta2.3 <- mean((par.est3[ , 2] - b2) ^2) # Mean squared error
perf.meas[3,4] <- mse.beta2.3

# Prediction error
pred.error3 <- mean(par.est3[ , 3])
perf.meas[3,5] <- pred.error3

# Make all these predictions into a nice dataframe
performance <- data.frame(perf.meas)
variablenames <- c("Absolute Bias b1", "Absolute Bias b2", "Mean Squared Error b1", "Mean Squared Error b2", "Prediction Error")
colnames(performance) <- variablenames
