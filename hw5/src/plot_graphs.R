# Author: Eric Kalosa-Kenyon
# Constructs and generates plots for the induced graphs in this assignment

## @knitr read_library

library(GGally)
library(network)

## @knitr graph1

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
plot.network(n, label=c("A", "B", "C", "D"),
             arrowhead.cex=2, label.cex=2, vertex.cex=2,
             vertex.col="green")
title("Figure 1 - from Bryan and Leise (2006)")

## @knitr graph2

add.vertices(n, 1)
add.edge(n, 5, 3)
add.edge(n, 3, 5)

plot.network(n, label=c("A", "B", "C", "D", "E"),
             arrowhead.cex=2, label.cex=2, vertex.cex=2,
             vertex.col="blue")
title("Figure 2 - Adding a new page")

## @knitr graph3

a2 = matrix(0, 4, 4)
a2[1,2] = 1;
a2[2,1] = 1;
a2[2,3] = 1;
a2[3,2] = 1;
a2[3,4] = 1;
a2[4,3] = 1;
n2 = network(a2, directed=T)
plot.network(n2, label=c("A", "B", "C", "D"),
             arrowhead.cex=2, label.cex=2, vertex.cex=2,
             vertex.col="red")
title("Figure 3 - a chain graph")
