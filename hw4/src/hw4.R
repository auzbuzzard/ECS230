# Author: Eric Kalosa-Kenyon
#
# Analysis script for polynomial fit

## @knitr READ
library(ggplot2)
library(dplyr)
library(glmnet)

# subroutines
html_table_width <- function(kable_output, width){
  width_html <- paste0(paste0('<col width="', width, '">'), collapse = "\n")
  sub("<table>", paste0("<table>\n", width_html), kable_output)
}

# parameters
data_dir = "../data/"

d = data.frame(D = numeric(), X = numeric(), Y = numeric())
# D == 0 is raw observations, otherwise D is degree of polynomial
for(i in 2:9){
    tmp = read.table(paste(data_dir, "poly_fit_", i, ".dat", sep=""))
    names(tmp) = c("X", "Y")
    tmp$D = i
    d = rbind(d, tmp)
}
tmp = read.table(paste(data_dir, "data.dat", sep=""), skip=1)
names(tmp) = c("X", "Y")
tmp$D = 0
d = rbind(d, tmp)
d$D = as.factor(d$D)

## @knitr GRAPH
plt0 = ggplot(d[d$D == 0,], aes(X, Y)) + geom_point() +
    labs(x="X", y="Y", title="Observed data")
plt0
plt1 = ggplot(d[d$D != 0,], aes(X, Y, color=D)) +
    geom_line() + facet_wrap(~ D) +
    labs(x="X", y="Y", title="Polynomial regressions of differing degrees",
         color = "Degree\n")
plt1

## @knitr POLY59
p5 = read.table("../data/poly_coef_5.dat")
names(p5) = c("degree.5.coefficients")
p9 = read.table("../data/poly_coef_9.dat")
names(p9) = c("degree.9.coefficients")

## @knitr AIC
dd = d[d$D == 0, c("X", "Y")]
a = data.frame(degree=numeric(), aic=numeric())
for(i in 1:8){
    m = lm(Y ~ poly(X, i), data=dd)
    # print(paste("AIC for poly of degree", d, ":", AIC(m)))
    a[i,] = c(i-1, AIC(m))
}

## @knitr LASSO
x = matrix(model.matrix(~ poly(X, 9), dd), 10, 10)[,2:10]
y = dd$Y
mlas = glmnet(x, y, alpha=1) # alpha=1 for LASSO
plot(mlas, xvar="lambda")

# crossval to select best lambda
cv.mlas = cv.glmnet(x, y, alpha=1) # alpha=1 for LASSO
plot(cv.mlas)

best.lambda <- cv.mlas$lambda.min

mlas.opt = glmnet(x, y, alpha=best.lambda)
