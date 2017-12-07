## @knitr graph1

library(GGally)
library(network)
a = matrix(0, 4, 4)
a[1,2] = 1;
a[1,3] = 1;
a[1,4] = 1;
a[2,3] = 1;
a[2,4] = 1;
a[3,1] = 1;
a[4,1] = 1;
a[4,3] = 1;
n = network(a, directed=T)
plot(n)
# ggnet(n, arrow.gap = 0.05, arrow.size=10)
