# http://blog.phytools.org/2018/01/adding-geological-legend-to-fan-style.html

# phytools type plot

tree$root.time <- max(ape::branching.times(tree))*1.2
strap::geoscalePhylo(tree = ape::ladderize(tree,right=FALSE),
                     units=c("Period"),
                     boxes="Period",
                     cex.tip=0.5,
                     cex.age=0.7,
                     cex.ts=0.7,
                     label.offset=0,
                     lwd=3,
                     width=2)

geoscalePhylo(tree=tree, boxes="Age", cex.tip=0.4)
