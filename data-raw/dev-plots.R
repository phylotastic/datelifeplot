# Adding packages to description
usethis::use_dev_package("phyloch", type = "imports", remote = "fmichonneau/phyloch")

# http://blog.phytools.org/2018/01/adding-geological-legend-to-fan-style.html

# ^ phytools type plot

# strap plot:
devtools::install("../datelifeplot/.")

phylo_sdm <- datelife::datelife_search(input = c("Delphinus_delphis", "Gallus gallus", "elephas Maximus", "felis_catus", "homo-sapiens"),
                                       use_tnrs = TRUE,
                                       summary_format = "phylo_sdm")
class(phylo_sdm) <- "phylo"

tree <- phylo_sdm
phylo_length <- max(ape::branching.times(tree))
time_depth <- round(phylo_length*1.2, digits = -1)
tree$root.time <- phylo_length
unit = "Period"
strap::geoscalePhylo(tree = tree,
                     x.lim = c(0, phylo_length),
                     cex.tip = 0.7,
                     show.tip.label = FALSE,
                     cex.ts = 0.7,
                     cex.age = 0.7,
                     width = 4,
                     tick.scale = "no",
# creates boxes with the last unit in argument "unit":
                     boxes = unit[length(unit)],
                     quat.rm = TRUE,
                     units = unit)


#
datelifeplot::plot_phylo(phylo_sdm, title = "", plot_type = "strap")
