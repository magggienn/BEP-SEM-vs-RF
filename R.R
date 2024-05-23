# 
setwd("/Users/magggien/Documents/BEP-SEM-vs-RF")
library(mvtnorm)
library(lavaan)
library(glmnet)
library(ggplot2)
library(colorspace)
library(randomForest)
library(corrr)
library(ggcorrplot)
library(factoextra)
library(MultivariateRandomForest)

rm(list=ls())
source("predicty.lavaan.R")

#define simulation
simstudy = function(samplesize = c(100, 200, 500, 1000), 
    measurement = c("M:weak", "M:strong"),
    structural = c("S:weak", "S:strong"),
    n2 = 1000,
    #skewdness of the data distribution
    skewtheta = FALSE,
    #
    direct = FALSE){

    # define lavaan model
    model <- '
    # latent variable definitions
    l1 =~ x1 + x2 + x3 + x4
    l2 =~ x5 + x6 + x7 + x8
    l3 =~ x9 + x10 + x11 + x12
    l4 =~ x13 + x14 + x15 + x16
    l5 =~ x17 + x18 + x19 + x20
    l6 =~ x21 + x22 + x23 + x24
    l7 =~ x25 + x26 + x27 + x28
    l8 =~ x29 + x30 + x31 + x32
    l9 =~ x33 + x34 + x35 + x36
    y =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + y11 + y12 + y13 + y14 + y15  + y16
    # regressions
    y ~ l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9
    '
    
    # variance covariance matrix of latent variables from CERQ study
    vcov.theta = matrix(c(
    1.000, 0.467, 0.425, 0.153, 0.465, 0.455, 0.445, 0.283, 0.101,
    0.467, 1.000, 0.471, 0.356, 0.493, 0.492, 0.513, 0.355, 0.256,
    0.425, 0.471, 1.000, 0.109, 0.652, 0.418, 0.142, 0.531, 0.327,
    0.153, 0.356, 0.109, 1.000, 0.324, 0.416, 0.564, 0.022, 0.063,
    0.465, 0.493, 0.652, 0.324, 1.000, 0.811, 0.430, 0.104, 0.114,
    0.455, 0.492, 0.418, 0.416, 0.811, 1.000, 0.685, -0.145, -0.070,
    0.445, 0.513, 0.142, 0.564, 0.430, 0.685, 1.000, -0.070, 0.060,
    0.283, 0.355, 0.531, 0.022, 0.104, -0.145, -0.070, 1.000, 0.664,
    0.101, 0.256, 0.327, 0.063, 0.114, -0.070, 0.060, 0.664, 1.000
    ), 9, 9)
    
    signbetas = sign(matrix(c( 2.712,
    -0.136,
    0.048,
    0.513,
    0.525,
    -2.671,
    0.042,
    5.324,
    -0.461), 9, 1))
    
    repetitions = 1
    
    if(!direct){
    PE = data.frame(repetition = rep(1:repetitions, each = 32),
    N = rep(rep(c(100, 200, 500,1000), each = 8), repetitions),
    M = rep(rep(c("M:weak", "M:strong"), each = 4), 2 * repetitions), 
    S = rep(rep(c("S:weak", "S:strong"), each = 2), 4 * repetitions),
    model = rep(c("sem", "rf"), 16 * repetitions),
    pe = rep(0, 32 * repetitions), 
    vareta = rep(0, 32 * repetitions))
    }
    else{
    PE = data.frame(repetition = rep(1:repetitions, each = 32),
    N = rep(rep(c(100, 200, 500,1000), each = 8), repetitions),
    M = rep(rep(c("M:weak", "M:strong"), each = 4), 2 * repetitions), 
    S = rep(rep(c("D:weak", "D:strong"), each = 2), 4 * repetitions),
    model = rep(c("sem", "rf"), 16 * repetitions),
    pe = rep(0, 32 * repetitions), 
    vareta = rep(0, 32 * repetitions))
    }
    
    teller1 = 0                
    set.seed(1234)
    
    for(rep in 1:repetitions){
    for(ss in 1:4){
    for(m in 1:2){
    for(s in 1:2){
    cat("This is repetition:", rep, "from", repetitions, "\n")
    teller1 = teller1 + 1
    
    n1 = samplesize[ss]
    meas = measurement[m]
    struc = structural[s]
    
    n = n1 + n2
    idx1 = 1:n1
    idx2 = (n1 + 1): n
    
    foldid = sample(rep(1:10, length.out = n1)) # for glmnet
    
    # data generation
    vareta = 1.1
    
    while(vareta > 1){
    # latent variables
    theta = rmvnorm(n, rep(0,9), vcov.theta)
    if(skewtheta){
    theta = scale(exp(theta)) # skew distribution of latent variables
    }
    
    #indicators 
    X = matrix(NA, n, 36)
    Y = matrix(NA, n, 16)
    
    # can adjust factor loadings for strong/weak measurement model
    if(meas == "M:strong"){
    lambda = runif(36, min = 0.5, max = 0.8)
    }
    else if(meas == "M:weak"){
    lambda = runif(36, min = 0.2, max = 0.5)
    } 
    if(!direct){
    if(struc == "S:strong"){
    betas = signbetas * matrix(runif(9, min = 0.25, max = 0.40), 9, 1)
    }
    else if(struc == "S:weak"){
    betas = signbetas * matrix(runif(9, min = 0.15, max = 0.25), 9, 1)
    }
    }
    else{
    betas = signbetas * matrix(runif(9, min = 0.15, max = 0.25), 9, 1)
    }
    
    # means
    mx = runif(36, min = 1.50, max = 3.0)
    my = runif(16, min = 1.15, max = 2.15)
    
    j = 0
    for(f in 1:9){
    for(r in 1:4){
    j = j + 1
    X[ , j] = mx[j] + lambda[j] * theta[, f] + rnorm(n, mean = 0, sd = sqrt(1 - lambda[j]^2))
    }
    }
    
    if(!direct){
    eta = theta %*% betas
    }
    else{
    dbetas = matrix(0, ncol(X), 1)
    if(struc == "S:strong"){
    dbetas[seq(3, 36, by = 4), 1] = (runif(9, min = 0.15, max = 0.25) * sample(c(-1,1), 9, replace = TRUE))
    }
    else if(struc == "S:weak"){
    dbetas[seq(3, 36, by = 4), 1] = (runif(9, min = 0.10, max = 0.15) * sample(c(-1,1), 9, replace = TRUE))
    }
    eta = theta %*% betas + X %*% dbetas
    }
    vareta = var(eta)
    }
    
    y = eta + rnorm(n, mean = 0, sd = sqrt(1 - vareta))
    
    lambday = runif(16, min = 0.2, max = 0.7)
    
    for(r in 1:16){
    Y[ , r] = my[r] + lambday[r] * eta + rnorm(n, mean = 0, sd = sqrt(1 - lambday[r]^2))
    }
    
    #final dataset 
    df = data.frame(cbind(X, Y))
    colnames(df) = c(paste("x", c(1:36), sep = ""), paste("y", c(1:16), sep = ""))
    xnames = paste("x", c(1:36), sep = "")
    ynames = paste("y", c(1:16), sep = "")
    dfpredictors = data.frame(X)
    dfoutcomes = data.frame(Y)
    
    rmse_list <- list()
    all_predictions_rf <- list()
    
    # Prepare the training and testing datasets
    train_data_rf <- dfpredictors[idx1, ]
    test_data_rf <- dfpredictors [idx2,]
    train_labels_rf <- dfoutcomes[idx1,]
    test_labels_rf <- dfoutcomes[idx2, ]
    
    n_tree <- 500 # Number of trees
    mtree <- 1 # Number of variables to possibly split at in each node
    min_leaf <- 1 # Minimum number of observations in each leaf node
    
    predictions_rf <- build_forest_predict(train_data_rf, train_labels_rf, n_tree, mtree, min_leaf, test_data_rf)
    write.csv(predictions_rf, file = "predictions_rf.csv")

    # Initialize lists to store results
    multirmse_list <- list()
    all_predictions_rf <- list()

    # Fit the Multivariate Random Forest model
    # predictions_rf <- MultivariateRandomForest(trainX, trainY, n_tree, mtree, min_leaf, testX)

    # Function to calculate RMSE for each outcome variable
    calculate_rmse <- function(true_values, predicted_values) {
      sqrt(mean((true_values - predicted_values)^2))
    }

    # Calculate RMSE for each outcome variable and store the predictions
    for (i in 1:ncol(test_labels_rf)) {
      rmse <- calculate_rmse(test_labels_rf[, i], predictions_rf[, i])
      multirmse_list[[paste("y", i, sep = "")]] <- rmse
      all_predictions_rf[[paste("y", i, sep = "")]] <- predictions_rf[, i]
    }

    # Print RMSE for each outcome variable
    print(multirmse_list)
    # 
    # 
    # for (i in 1:ncol(dfoutcomes)){
    #   lv_column <- dfoutcomes[, i]
    #   df_rf = data.frame(cbind(dfpredictors, lv_column))
    # 
    #   
    #   # Prepare the training and testing datasets
    #   train_data_rf <- df_rf[idx1, -ncol(df_rf)]
    #   test_data_rf <- df_rf[idx2, -ncol(df_rf)]
    #   train_labels_rf <- df_rf[idx1, ncol(df_rf)]
    #   test_labels_rf <- df_rf[idx2, ncol(df_rf)]
    #   
    #   # Fit the Random Forest model
    #   rf_model <- randomForest(x = train_data_rf, y = train_labels_rf, ntree = 500)
    #   yhat_rf <- predict(rf_model, newdata = test_data_rf)
    #   
    #   # Add the predictions to the list with the column name
    #   all_predictions_rf[[colnames(dfoutcomes)[i]]] <- yhat_rf
    #   
    #   # Print the RMSE for each outcome variable
    #   rmse_rf <- sqrt(mean((test_labels_rf - yhat_rf)^2))
    #   cat("Outcome Variable", i, "RMSE:", rmse_rf, "\n")
    #   # Add the RMSE to the list with the column name as key
    #   rmse_list[[colnames(dfoutcomes)[i]]] <- rmse_rf
    # }
    # 
    # # Combine RMSEs into a single data frame
    # rmse_rf <- data.frame(variable = names(rmse_list), rmse = unlist(rmse_list), row.names = NULL)

  
    # fit SEM and predict
    fit <- sem(model, data = df[idx1, ], std.lv = TRUE, meanstructure = TRUE, warn = FALSE)
    yhat.sem = predicty.lavaan(fit, newdata = df[idx2, ], xnames = xnames, ynames = ynames)

    # store results
    pe.iter = c(sqrt(mean((Y[idx2,] - yhat.sem)^2)),multirmse_list$rmse)
    
    # pe.iter = c(sqrt(mean((Y[idx2,] - yhat.sem)^2)))
    PE$pe[((teller1 -1)*2 + 1): (teller1*2)] = pe.iter
    PE$vareta[((teller1 -1) * 2 + 1): (teller1 * 2)] = rep(vareta, 2)
    
    }
    }
    }
    }
    PE$M = as.factor(PE$M)
    PE$N = as.factor(PE$N)
    PE$S = as.factor(PE$S)
    PE$model = as.factor(PE$model)
    
    return(PE)
    
}

PE1 = simstudy(samplesize = c(100, 200, 500, 1000), 
                   measurement = c("M:weak", "M:strong"),
                   structural = c("S:weak", "S:strong"),
                   n2 = 1000,
                   skewtheta = FALSE, 
                   direct = FALSE)

PE2 = simstudy(samplesize = c(100, 200, 500, 1000), 
               measurement = c("M:weak", "M:strong"),
               structural = c("S:weak", "S:strong"),
               n2 = 1000,
               skewtheta = TRUE, 
               direct = FALSE)

PE3 = simstudy(samplesize = c(100, 200, 500, 1000), 
               measurement = c("M:weak", "M:strong"),
               structural = c("S:weak", "S:strong"),
               n2 = 1000,
               skewtheta = FALSE, 
               direct = TRUE)


# save results
save(PE1, file = "simstudy1final_newrf_16.Rdata")
save(PE2, file = "simstudy2final.Rdata")
save(PE3, file = "simstudy3final_newrf_16.Rdata")

# plot results
mycolor = c(sequential_hcl(5, palette = "Reds 3")[2], sequential_hcl(5, palette = "Greens 2")[2])

p1 <- ggplot(PE1, aes(x=model, y=pe, fill=factor(model))) + 
  geom_boxplot(aes(group = factor(model))) + 
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) + 
  facet_grid(N ~ M * S) +
  xlab("Model") + 
  ylab("RMSEp") + 
  scale_fill_manual(values = mycolor) + 
  ggtitle("Simulation Study 1") 
#theme_apa(legend.pos = "none")

ggsave("/Users/magggien/Documents/BEP-SEM-vs-RF/sim1_rf_10rep_newrf.pdf", 
       plot = p1, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
print(p1)

p2 <- ggplot(PE2, aes(x=model, y=pe, fill=factor(model))) + 
  geom_boxplot(aes(group = factor(model))) + 
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) + 
  facet_grid(N ~ M * S) +
  xlab("Model") + 
  ylab("RMSEp") + 
  scale_fill_manual(values = mycolor) + 
  ggtitle("Simulation Study 2") 
#+theme_apa(legend.pos = "none")
print(p2)

ggsave("/Users/magggien/Documents/BEP-SEM-vs-RF/sim2_test_workspacefile.pdf", 
       plot = p2, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)

p3 <- ggplot(PE3, aes(x=model, y=pe, fill=factor(model))) + 
  geom_boxplot(aes(group = factor(model))) + 
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) + 
  facet_grid(N ~ M * S) +
  xlab("Model") + 
  ylab("RMSEp") + 
  scale_fill_manual(values = mycolor) + 
  ggtitle("Simulation Study 3")  
#+theme_apa(legend.pos = "none")

ggsave("/Users/magggien/Documents/BEP-SEM-vs-RF/sim3_rfnew.pdf", 
       plot = p3, width = 11.7, height = 8.3, units = "in", limitsize = FALSE)
print(p3)

# make summaries
library(tidyverse)
PE1 %>%
  group_by(N, M, S) %>%
  summarize(mean = mean(vareta))

sum.pe1 = PE1 %>%
  group_by(N, M, S, model) %>%
  summarize(mean = mean(pe))

sum.pe2 = PE2 %>%
  group_by(N, M, S, model) %>%
  summarize(mean = mean(pe))

sum.pe3 = PE3 %>%
  group_by(N, M, S, model) %>%
  summarize(mean = mean(pe))

PE3 %>%
  group_by(N, M, S) %>%
  summarize(mean = mean(vareta))

# # ANOVAs
# library(rstatix)
# library(xtable)
# 
# load("/Users/magggien/Documents/BEP-SEM-vs-RF/simstudy1final.Rdata")
# PE1$id = rep(1:16, each = 2)
# res.aov1 <- anova_test(
#   data = PE1, 
#   dv = pe, 
#   wid = id,
#   within = model, 
#   between = c(N, M, S),
#   effect.size = "pes")
# xtable(res.aov1)
# 
# load("/Users/magggien/Documents/BEP-SEM-vs-RF/simstudy2final_newrf_16.Rdata")
# PE2$id = rep(1:16, each = 2)
# res.aov2 <- anova_test(
#   data = PE2, 
#   dv = pe, 
#   wid = id,
#   within = model, 
#   between = c(N, M, S),
#   effect.size = "pes")
# xtable(res.aov2)
# 
load("/Users/magggien/Documents/BEP-SEM-vs-RF/simstudy1final_newrf_16.Rdata")
PE3$id = rep(1:16, each = 2)
res.aov3 <- anova_test(
  data = PE3,
  dv = pe,
  wid = id,
  within = model,
  between = c(N, M, S),
  effect.size = "pes")
xtable(res.aov3)

