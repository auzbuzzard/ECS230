# Author: Eric Kalosa-Kenyon
#
# Analysis script for repeated matrix multiplications

## @knitr READ
df = read.table("n50_r100_ord1.csv", header=T, sep=",")
head(df)

## @knitr PLOT
plot(df)
hist(df$clocks)
hist(df$time)
