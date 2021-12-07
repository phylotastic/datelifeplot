# Adding packages to description
usethis::use_dev_package("phyloch", type = "imports", remote = "fmichonneau/phyloch")

# http://blog.phytools.org/2018/01/adding-geological-legend-to-fan-style.html

# ^ phytools type plot

# strap plot:
devtools::install("../datelifeplot/.")
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

datelifeplot::plot_phylo(phylo_sdm, title = "", axis_type = 3)

strap::geoscalePhylo(tree = tree,
                     x.lim = c(0, time_depth),
                     cex.tip = 0.7,
                     cex.ts = 0.7,
                     cex.age = 0.7,
                     width = 4,
                     tick.scale = 15,
                     boxes = "Epoch",
                     erotate = 90,
                     quat.rm = TRUE,
                     units = c("Period","Epoch"))
