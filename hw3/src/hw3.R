# Author: Eric Kalosa-Kenyon
#
# Analysis script for repeated matrix multiplications

## @knitr SETUP

data_dir = "../data/"
simfiles_ord1 = c("n100_r1000_ord1.csv", "n200_r100_ord1.csv",
                  "n500_r20_ord1.csv", "n1000_r20_ord1.csv",
                  "n2000_r10_ord1.csv")

simfiles_loop = c("n250_r50_ordAll.csv")

simfiles_opt2 = c("n100_r1000_ordAll_O2opt.csv", "n250_r500_ordAll_O2opt.csv",
                  "n500_r200_ordAll_O2opt.csv", "n2000_r100_ordAll_O2opt.csv")

simfiles_opt3 = c("n100_r1000_ordAll_O3opt.csv", "n250_r500_ordAll_O3opt.csv",
                  "n500_r200_ordAll_O3opt.csv", "n2000_r100_ordAll_O3opt.csv")

## @knitr READ_pr1
df1 = data.frame(clocks=NA, time=NA, batch=NA)
for(fn in simfiles_ord1){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$batch = fn
    df1 = rbind(df1, df)
}
df1 = na.omit(df1)

## @knitr READ_pr2
simfiles_pr2 = c(simfiles_loop, simfiles_opt2, simfiles_opt3)
df2 = data.frame(clocks=NA, time=NA, loop_order=NA, batch=NA)
for(fn in simfiles_pr2){
    fp = paste(data_dir, fn, sep="")
    df = read.table(fp, header=T, sep=",")
    df$batch = fn
    df2 = rbind(df2, df)
}
df2 = na.omit(df2)


## @knitr PLOT_pr1
# plot(df)
# hist(df$clocks)
# hist(df$time)
