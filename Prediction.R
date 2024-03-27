fit = sem(mymodel, data = train, meanstructure = TRUE)
yhat = preditcty.lavaan (fit, newdata = test, xnames = xn, ynames = yn) #xnames - predictor variables, ynames - response varbia