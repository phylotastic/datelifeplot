# the node position on the y axis is given by pp$yy
# pp$yy also stores the positions of the tips
# the y coordinate of nodes start at pp$yy[ape::Ntip + 1]
#' Add Node Age/Calibration Data to a Plotted Phylogeny
#' @param phy A phylo object with branch length proportional to time
#' @param cex_tips A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of `phy` tip labels. Default to 1.
#' @param calibration_summary An output of [datelife:::summary.matchedCalibrations()]
#' @param color_points A named vector of colors. Names must correspond to study
#'  in `calibration_summary$references`.
#'  If vector is not named, colors will be recycled.
#' @param pch A numeric vector indicating the type of point to plot.
#'   See [graphics::par()] for options. Defaults to 20 = "bullet circle".
#' @param cex_pch A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of node age point symbols defined by `pch`. Default to 1.
#' @param color_bars A character vector of one element indicating the color for
#'  the bars. Default to "gray".
#' @param lty_bars A numeric vector or character string indicating the plotting
#'   line type for bars. See [graphics::par()] for options. Defaults to 1 = "solid".
#' @param lwd_bars A numeric vector indicating the plotting width for bars.
#'   See [graphics::par()] for options. Defaults to 2.
#' @importFrom ape .PlotPhyloEnv
#' @export
plot_node_ages <- function(phy,
                           cex_tips = 1,
                           calibration_summary,
                           color_points,
                           pch = 20, # bullet circle
                           cex_pch = 1,
                           color_bars,
                           lty_bars,
                           lwd_bars) {

  # obtain max x lim from ages
  x_max <- max(c(calibration_summary$in_phy$MinAge, calibration_summary$in_phy$MaxAge))
  # plot.phylo
  ape::plot.phylo(phy,
                  cex = cex_tips,
                  plot = FALSE)
  # get plotting x.lim
  lastPP_x.lim <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2]
  ape::plot.phylo(phy,
                  cex = cex_tips,
                  x.lim = c(-5,lastPP_x.lim))
  # get recorded plotting parameters
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  # we use lastPP$xx positions to get plot x position for node ages
  if (all(calibration_summary$in_phy$MinAge == calibration_summary$in_phy$MaxAge)) {
    # We need the one max number in lastPP$xx, then we substract the node ages to get
    # the actual x positions for all node ages on the plot:
    x_ages <- max(lastPP$xx) - calibration_summary$in_phy$MinAge
  }
  # we use lastPP$yy positions to get plotting y position of nodes
  if (all(calibration_summary$in_phy$mrca_node_number > lastPP$Ntip)) {
    # case in which mrca_node_numbers start at Ntip + 1
    # we can directly use mrca_node_numbers as index for lastPP$yy
    # it generates a vector of the correct length
    y_nodes <- lastPP$yy[calibration_summary$in_phy$mrca_node_number]
  }
  # using references to color points
  if (missing(color_points)) {
    color_points <- calibration_summary$in_phy$reference
  }
  for (i in unique(calibration_summary$in_phy$mrca_node_number)) {
    rowsies <- calibration_summary$in_phy$mrca_node_number == i
    x_min <- min(x_ages[rowsies])
    x_max <- max(x_ages[rowsies])
    # we just need one element for y position of nodes:
    y_pos <- y_nodes[rowsies][1]
    graphics::segments(x0 = x_min,
             y0 = y_pos,
             x1 = x_max,
             y1 = y_pos,
             col = ifelse(missing(color_bars), "gray", color_bars),
             lty = ifelse(missing(lty_bars), "solid", lty_bars),
             lwd = ifelse(missing(lwd_bars), 2, lwd_bars))
  }
  graphics::points(x_ages,
         y_nodes,
         col = color_points,
         pch = pch)
  ape::axisPhylo(1)
  graphics::mtext("Time (myrs)", side = 1, line = 2, at = max(lastPP$xx) * 0.5)
}
