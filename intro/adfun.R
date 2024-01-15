library(RTMB)

f<-function(x){
  x[1]^2+sqrt(x[2])+x[1]^2*sqrt(x[2])
}

F<-MakeTape(f, c(2,3))

F

F(c(1,2))

F$jacobian(c(1,2))

DF<-F$jacfun()

DF

DF$jacobian(c(1,2))

F$graph()

library(igraph)
G <- graph_from_adjacency_matrix(F$graph())
plot(G, vertex.size=17, layout=layout_as_tree)
