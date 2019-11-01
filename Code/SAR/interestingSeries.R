x = c(1:100)
qplot(x = c(1:100), y = sin(2*x)*cos(2*x), geom = "line", main = expression(paste("sin(", 2 * x, ")", "cos(", 2 * x, ")")), xlab = "X", ylab = "Y") +
  theme_few(base_size = 18, base_family = "serif")  + theme(plot.title = element_text(hjust=0.5)) + 
  scale_colour_few("Dark")

ts =  sin(2*x)*cos(2*x)
ts = ts/max(ts)
e = sliding.window(ts, 3, 1)
p = formation.pattern(e)
w = pattern.wedding(p)
g = transition.graph(e, w, 3)

m = as.matrix(g)
net=graph.adjacency(adjmatrix=m,mode="directed",weighted=TRUE,diag=TRUE)


V(net)$color <- "#C9C5C5"  
V(net)$size <- 40
E(net)$arrow.size <- .1
E(net)$edge.color <- "#838181"
E(net)$width <- 1
E(net)$label <- round(E(net)$weight, 2)
plot(net, vertex.label = V(net)$name)
