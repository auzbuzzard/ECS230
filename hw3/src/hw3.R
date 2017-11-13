# Author: Eric Kalosa-Kenyon
#
# Analysis script for repeated matrix multiplications

## @knitr SETUP
library(ggplot2)

data_dir = "../data/"

simfiles_ord1 = c("n100_r1000_ord1.csv", "n200_r100_ord1.csv",
                  "n500_r20_ord1.csv", "n1000_r20_ord1.csv",
                  "n2000_r10_ord1.csv")

simfiles_loop = c("n250_r50_ordAll.csv")

simfiles_opt2 = c("n100_r1000_ordAll_O2opt.csv", "n250_r500_ordAll_O2opt.csv",
                  "n500_r200_ordAll_O2opt.csv", "n2000_r100_ordAll_O2opt.csv")

simfiles_opt3 = c("n100_r1000_ordAll_O3opt.csv", "n250_r500_ordAll_O3opt.csv",
                  "n500_r200_ordAll_O3opt.csv", "n2000_r100_ordAll_O3opt.csv")

## @knitr PR1
# Read simulated matrix multiplications
df1 = data.frame(clocks=NA, time=NA, batch=NA)
for(fn in simfiles_ord1){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$batch = fn
    # remove outliers
    rmv.ixs = which(df$clocks > mean(df$clocks) + 2*sqrt(var(df$clocks)))
    df = df[-rmv.ixs,]
    df1 = rbind(df1, df)
}
df1 = na.omit(df1)

# Plot matrix multiplication timing
plt1 = ggplot(df1, aes(clocks, fill=batch)) + geom_histogram() +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
plt1
plt2 = ggplot(df1, aes(time, fill=batch)) + geom_histogram() +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
plt2

## @knitr PR2
# Read data
simfiles_pr2 = c(simfiles_loop, simfiles_opt2, simfiles_opt3)
df2 = data.frame(clocks=NA, time=NA, loop_order=NA, batch=NA)
for(fn in simfiles_pr2){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$batch = fn
    # remove outliers
    rmv.ixs = which(df$clocks > mean(df$clocks) + 2*sqrt(var(df$clocks)))
    df = df[-rmv.ixs,]
    df2 = rbind(df2, df)
}
df2 = na.omit(df2)
df2$loop_order = factor(df2$loop_order)

# Plot matrix multiplication timing
plt3 = ggplot(df2[df2$batch == df2$batch[1],], # use only one batch at first
              aes(clocks, fill=loop_order)) + geom_histogram() +
    facet_wrap(~loop_order, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
plt3
plt4 = ggplot(df2, aes(clocks, fill=loop_order)) + geom_histogram() +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
plt4
