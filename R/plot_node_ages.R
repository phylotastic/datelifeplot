# the node position on the y axis is given by pp$yy
# pp$yy also stores the positions of the tips
# the y coordinate of nodes start at pp$yy[ape::Ntip + 1]
#' Add Node Age/Calibration Data to a Plotted Phylogeny
#' @param phy A phylo object with branch length proportional to time
#' @param cex_tips A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of `phy` tip labels. Default to 1.
#' @param calibration_summary An output of [datelife:::summary.matchedCalibrations()]
#' @param color A named vector of colors. Names must correspond to study references.
#'       If vector is not named, colors will be recycled.
#' @param cex_pch A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of node age point symbols defined by `pch`. Default to 1.
#' @export
plot_node_ages <- function(phy,
                           cex_tips = 1,
                           calibration_summary,
                           color,
                           pch = 20, # bullet circle
                           cex_pch = 1) {

  # obtain max x lim from ages
  x_max <- max(c(calibration_summary$in_phy$MinAge, calibration_summary$in_phy$MaxAge))
  # plot.phylo
  ape::plot.phylo(phy,
                  cex = cex_tips,
                  plot = FALSE)
  # get plotting x.lim
  library(ape) # we need this so it loads object .PlotPhyloEnv
  lastPP_x.lim <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2]
  ape::plot.phylo(phy,
                  cex = cex_tips,
                  x.lim = c(-5,lastPP_x.lim))
  # get recorded plotting parameters
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  # we use lastPP$xx positions to get plotting coordinates for node ages
  if (all(calibration_summary$in_phy$MinAge == calibration_summary$in_phy$MaxAge)) {
    x_ages <- max(lastPP$xx) - calibration_summary$in_phy$MinAge
  }
  # we use lastPP$yy positions to get plotting coordinates for node numbers
  if (all(calibration_summary$in_phy$mrca_node_number > lastPP$Ntip)) {
    # case in which mrca_node_numbers start at Ntip + 1
    # we can directly use mrca_node_numbers as index for lastPP$yy
    y_nodes <- lastPP$yy[calibration_summary$in_phy$mrca_node_number]
  }
  # using references to color points
  if (missing(color)) {
    color <- calibration_summary$in_phy$reference
  }
  points(x_ages,
         y_nodes,
         col = color,
         pch = pch)
  ape::axisPhylo(1)
  graphics::mtext("Time (myrs)", side = 1, line = 2, at = max(lastPP$xx) * 0.5)
}
