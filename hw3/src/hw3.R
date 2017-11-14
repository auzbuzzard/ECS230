# Author: Eric Kalosa-Kenyon
#
# Analysis script for repeated matrix multiplications

## @knitr SETUP
library(ggplot2)

data_dir = "../data/"

simfiles_ord1 = c("n1000_r20_ord1.csv", "n100_r1000_ord1.csv",
                  "n2000_r10_ord1.csv", "n200_r100_ord1.csv",
                  "n500_r20_ord1.csv")

simfiles_loop = c("n100_r100_ordAll.csv", "n250_r100_ordAll.csv")

simfiles_opt2 = c("n100_r1000_ordAll_O2opt.csv", "n250_r500_ordAll_O2opt.csv",
                  "n500_r200_ordAll_O2opt.csv", "n2000_r100_ordAll_O2opt.csv")

simfiles_opt3 = c("n1000_r50_ordAll_O3opt.csv", "n100_r1000_ordAll_O3opt.csv",
                  "n2000_r20_ordAll_O3opt.csv", "n250_r250_ordAll_O3opt.csv",
                  "n500_r100_ordAll_O3opt.csv")

simfiles_dgemm = c("n1000_r100_dgemm.csv", "n100_r1000_dgemm.csv",
                   "n2000_r50_dgemm.csv", "n200_r500_dgemm.csv",
                   "n500_r250_dgemm.csv")

get_batch = function(fnm){
    r = strsplit(strsplit(fnm, "_")[[1]][1], "n")[[1]][2]
    return(r)
}

## @knitr PR1
# Read simulated matrix multiplications
df1 = data.frame(clocks=NA, time=NA, fn=NA)
for(fn in simfiles_ord1){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$fn = fn
    # remove outliers
    # rmv.ixs = which(df$clocks > mean(df$clocks) + 4*sqrt(var(df$clocks)))
    # df = df[-rmv.ixs,]
    df1 = rbind(df1, df)
}
df1 = na.omit(df1)
df1$batch = unlist(lapply(df1$fn, get_batch))
df1$batch = factor(df1$batch)

# Plot matrix multiplication timing
plt1 = ggplot(df1, aes(clocks, fill=batch)) +
    geom_histogram(bins=50) + geom_freqpoly(bins=50) +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    ggtitle("Timing matrix multiplication with different size matrices") +
    ylab("Frequency") + xlab("Clock cycles used") +
    scale_colour_discrete(name = "Matrix size (n by n)")
plt1

plt2 = ggplot(df1, aes(time, fill=batch)) +
    geom_histogram(bins=50) + geom_freqpoly(bins=50) +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    ggtitle("Timing matrix multiplication with different size matrices") +
    ylab("Frequency") + xlab("Seconds used") +
    scale_colour_discrete(name = "Matrix size (n by n)")
plt2

## @knitr PR2
# Read data
simfiles_pr2 = c(simfiles_loop) # , simfiles_opt2) #, simfiles_opt3)
df2 = data.frame(clocks=NA, time=NA, loop_order=NA, fn=NA)
for(fn in simfiles_pr2){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$fn = fn
    df2 = rbind(df2, df)
}
df2 = na.omit(df2)
df2$loop_order = factor(df2$loop_order)
df2$batch = unlist(lapply(df2$fn, get_batch))
df2$batch = factor(df2$batch)

# Plot matrix multiplication timing
plt3 = ggplot(df2[df2$batch == df2$batch[1],], # use only one batch at first
              aes(clocks, fill=loop_order)) +
    geom_histogram(bins=50) + geom_freqpoly(bins=50) +
    facet_wrap(~loop_order, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    ggtitle("Timing matrix multiplication with different loop orderings") +
    ylab("Frequency") + xlab("Clock cycles used")
plt3

plt4 = ggplot(df2, aes(clocks, fill=loop_order)) +
    geom_histogram(bins=50) +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    ggtitle("Timing matrix multiplication with different loop orderings") +
    ylab("Frequency") + xlab("Clock cycles used")
plt4

## @knitr PR3
df3 = data.frame(clocks=NA, time=NA, fn=NA)
for(fn in simfiles_dgemm){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$fn = fn
    df3 = rbind(df3, df)
}
df3 = na.omit(df3)
df3$batch = unlist(lapply(df3$fn, get_batch))
df3$batch = factor(df3$batch)

plt5 = ggplot(df3, aes(clocks, fill=batch)) +
    geom_histogram(bins=50) + geom_freqpoly(bins=50) +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none") +
    ggtitle("Timing matrix multiplication using BLAS dgemm()") +
    ylab("Frequency") + xlab("Clock cycles used")
plt5

plt6 = ggplot(df3, aes(time, fill=batch)) +
    geom_histogram(bins=50) + geom_freqpoly(bins=50) +
    facet_wrap(~batch, scales="free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position="none") +
    ggtitle("Timing matrix multiplication using BLAS dgemm()") +
    ylab("Frequency") + xlab("Seconds used")
plt6
