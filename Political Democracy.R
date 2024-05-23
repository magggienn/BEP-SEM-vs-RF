library(lavaan)
data(PoliticalDemocracy)
colnames(PoliticalDemocracy) = c("z1", "z2", "z3", "z4", "y1", "y2", "y3", "y4", "x1", "x2", "x3")
head(PoliticalDemocracy)

model0 <- '
# latent variable definitions
ind60 =~ x1 + x2 + x3
dem60 =~ z1 + z2 + z3 + z4
dem65 =~ y1 + y2 + y3 + y4
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
z1 ~~ y1
z2 ~~ z4 + y2
z3 ~~ y3
z4 ~~ y4
y2 ~~ y4
'
fit <- sem(model0, data = PoliticalDemocracy, meanstructure = TRUE, warn = FALSE)
semPaths(fit, title = FALSE, intercepts = FALSE, residuals = FALSE)

xnames = colnames(PoliticalDemocracy)[-c(5,6,7,8)]
ynames = colnames(PoliticalDemocracy)[c(5,6,7,8)]
set.seed(1234)
repeats = 100
PE = data.frame(repetition = rep(1:repeats, each = 2),
                model = rep(1:2, repeats),
                pe = rep(0, 2 * repeats))
folds = rep(1:10, length.out = 75)
t = 0
for (r in 1:repeats){
  yhat1 = yhat2 = matrix(NA, 75, 4)
  folds = sample(folds)
  for(k in 1:10){
    t = t + 1
    idx = which(folds == k)
    # SEM approach
    fit <- sem(model0, data = PoliticalDemocracy[-idx, ], meanstructure = TRUE, warn = FALSE)
    yhat1[idx, ] = predicty.lavaan(fit, newdata = PoliticalDemocracy[idx, ], xnames = xnames, ynames = ynames)
  }
  pe1 = sqrt(sum((PoliticalDemocracy[, ynames] - yhat1)**2)/300)
  PE$pe = c(pe1)
}


library(ggplot2)
PE$model = as.factor(PE$model)
p <- ggplot(PE, aes(x=model, y=pe, fill=factor(model))) +
  geom_boxplot(aes(group = factor(model))) +
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) +
  xlab("Approach") + ylab("RMSEp") +
  theme(legend.position="none") +
  scale_fill_grey(start=.3,end=.7)
p
