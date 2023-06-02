library(FossilSim)
set.seed(1234567888)
# beta=0 - speciation occurs though budding

#simulate complete tree
t = TreeSim::sim.bd.taxa(n=10, numbsim=1, lambda=0.5, mu=0., frac=1.)[[1]]
# simulate taxonomy using fossilsim
beta = 0. # probability of symmetric speciation
s = sim.taxonomy(tree = t, beta = beta) 
# simulate fossils
f = sim.fossils.poisson(rate=0.5, taxonomy=s)
#f = sim.extant.samples(f, taxonomy=s, rho = 0.5)
plot(f, t, taxonomy = s, show.ranges = T, show.taxonomy = T)
plot(t)

# transform format
t2 = SAtree.from.fossils(t,f)$tree
# transform to sampled tree
t3 = sampled.tree.from.combined(t2, rho = 1.)
plot(t3)


# simulate fossils
f = sim.fossils.poisson(rate=0.1, tree=t, taxonomy=)
# simulate extant samples
f = sim.extant.samples(f, t, rho = 0.5)
plot(f, t)

f = sim.fossils.poisson(0.1, t, root.edge = FALSE)
 
# plot the output
plot(s, tree = t, legend.position = "bottomright")

f = reconcile.fossils.taxonomy(f, s)
plot(f, t, taxonomy = s, show.ranges = T, show.taxonomy = T)


b=fossils.to.BEAST.start.tree(t, f, complete = T)
fossils.to.BEAST.constraints(f,t, complete=T)
beast.fbd.format(t,f)
# transform format
t2 = SAtree.from.fossils(t,f)$tree
# transform to sampled tree
t3 = sampled.tree.from.combined(t2)
plot(t3)

t<- sim.fbd.taxa(5,1,0.1,0.05,0.1,0.5)[[1]]
plot(t)
rangeplot.asymmetric(t)
