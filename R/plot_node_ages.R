# @details The node position on the y axis is given by pp$yy
# pp$yy also stores the positions of the tips
# the y coordinate of nodes start at pp$yy[ape::Ntip + 1]

#' Plot a Chronogram and Add Age Data to Nodes
#' @description `plot_node_ages` plots a chronograms given in `phy` and adds age
#'   data [points()] given in `calibration_summary` to the corresponding nodes.
#'
#' @param phy A `phylo` object with branch length proportional to time.
#' @param calibration_summary An output of [datelife:::summary.matchedCalibrations()]
#' @param color_pch A named vector of colors. Names must correspond to study names
#'  in `calibration_summary$in_phy$references`.
#'  If vector is not named, colors will be recycled.
#' @param transparent_pch A numeric value ranging from 10-99, indicating
#' transparency for `color_pch`. Default to `NULL`, no transparency.
#' @param pch A numeric vector indicating the symbol to plot age points.
#'   See [graphics::par()] for options. Defaults to 20 = "bullet circle".
#' @param cex_pch A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of node age point symbols defined by `pch`. Default to 1.
#' @param color_bars A character vector of one element indicating the color for
#'  node age distribution bars. Defaults to "#80808050", which is hex for gray
#' ("#808080") with 50\% transparency.
#' @param lty_bars A numeric vector or character string indicating the plotting
#'   line type for bars. See [graphics::par()] for options. Default to "solid".
#' @param lwd_bars A numeric vector indicating the line width for age distribution bars.
#'   See [graphics::par()] for options. Default to 7.
#' @param add_legend Default to `TRUE`, adds a legend to the left of the plot.
#' @param cex_legend A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of legend. Default to one half the size of the axis label, `cex_axislabel * 0.5`.
#' @param x_legend the x co-ordinate to be used to position the legend on the left side of the plot.
#' @param y_legend the y co-ordinate to be used to position the legend on the left side of the plot.
#' @param mai1,mai2,mai3,mai4 A numeric value indicating internal plot margin sizes ininches.
#' @param omi1,omi2,omi3,omi4 A numeric value indicating outter plot margin sizes in inches.
#' @inheritParams plot_phylo
#' @inheritDotParams ape::plot.phylo
#' @importFrom ape .PlotPhyloEnv
#'@details Plot are margin sizes as defined by [graphics::par()$mai] and [graphics::par()$omi]
#' are overruled within the function. To modify them you have to use the arguments
#' `mai1`, `mai2`, `mai3` and `mai4`, and omi1, omi2, omi3 and omi4.
#' @export
plot_node_ages <- function(phy,
                           #plotting_method = "plot_phylo",
                           time_depth = NULL,
                           plot_type = "phyloch",
                           mai1, mai2, mai3, mai4,
                           omi1, omi2, omi3, omi4,
                           title = "Chronogram",
                           cex_title = graphics::par("cex"),
                           pos_title = 1,
                           cex_tiplabels = graphics::par("cex"),
                           geologic_timescale = "strat2012",
                           geologic_unit = "period",
                           pos_axis = 1,
                           cex_axis = graphics::par("cex"),
                           axis_label = "Time (MYA)",
                           center_axislabel = 0.5,
                           cex_axislabel = graphics::par("cex"),
                           calibration_summary,
                           pch = 20,
                           color_pch,
                           transparent_pch = NULL,
                           cex_pch = graphics::par("cex"),
                           color_bars = "#80808050",
                           lty_bars = "solid",
                           lwd_bars = 7,
                           add_legend = TRUE,
                           cex_legend = cex_axislabel*0.5,
                           x_legend = NULL,
                           y_legend = NULL,
                           ...) {
  # get calibrations that are in_phy only
  in_phy <- calibration_summary$in_phy
  # obtain max x lim from calibration ages, chronogram and chronogram + root:
  phylo_length <- max(ape::branching.times(phy))
  max_calibration  <- max(c(in_phy$MinAge, in_phy$MaxAge))
  if (is.null(time_depth)) {
    print(paste("time_depth = ", time_depth))
    if (is.null(phy$root.edge)) {
      time_depth <- round(max(phylo_length, max_calibration)*1.2, digits = -1)
    } else {
      time_depth <- round(max(phylo_length + phy$root.edge, max_calibration), digits = -1)
    }
  }
  # add a root to the chronogram (or extend the root) to plot from the specified time_depth
  phy$root.edge <- time_depth - max(ape::branching.times(phy))
  # define plotting area margins
  if (missing(mai1)) {
    mai1 <- 0
  }
  if (missing(mai2)) {
    mai2 <- ifelse(add_legend, 2, 0)
  }
  if (missing(mai3)) {
    mai3 <- 0
  }
  if (missing(mai4)) {
    ind <- which.max(nchar(phy$tip.label))
    mai4 <- graphics::strwidth(s = phy$tip.label[ind],
                               units = "inches",
                               cex = cex_tiplabels,
                               font = 3)
  }
  if (missing(omi1)) {
    pho <- phylo_height_omi(phy = phy)
    omi1 <- pho$omi1
  }
  if (missing(omi2)) {
    omi2 <- 0
  }
  if (missing(omi3)) {
    omi3 <- 1
  }
  if (missing(omi4)) {
    omi4 <- 0
  }
  graphics::par(xpd = NA,
                mai = c(0, mai2, 0, mai4),
                # par()$mai[2]
                omi = c(omi1, 0, 1, 0))

  # plot.phylo
  ape::plot.phylo(phy,
                  cex = cex_tiplabels,
                  label.offset = 0.5,
                  x.lim = c(0, time_depth),
                  root.edge = TRUE,
                  plot = FALSE,
                  ...)
  # get plotting x.lim
  lastPP_x.lim <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2]
  ape::plot.phylo(phy,
                  cex = cex_tiplabels, #edge.width = 2,
                  label.offset = 0.5,
                  x.lim = c(0, lastPP_x.lim),
                  root.edge = TRUE,
                  plot = TRUE, ...)

  # get recorded plotting parameters
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  # we use lastPP$xx positions to get plot x position for node ages
  ages_equal <- in_phy$MinAge == in_phy$MaxAge
  if (all(ages_equal)) {
    # We need the one max x position in lastPP$xx with max(lastPP$xx)
    # Then we substract the node ages to get actual x positions for
    # all node ages on the plot:
    x_ages <- max(lastPP$xx) - in_phy$MinAge
  }
  # Next, we use lastPP$yy positions to get plotting y position of nodes
  if (all(in_phy$mrca_node_number > lastPP$Ntip)) {
    # case in which mrca_node_numbers start at Ntip + 1
    # we can directly use mrca_node_numbers as index for lastPP$yy
    # it generates a vector of the correct length
    y_nodes <- lastPP$yy[in_phy$mrca_node_number]
  }
  # Use study references to color points
  if (missing(color_pch)) {
    color_pch_all <- in_phy$reference
  } else {
    color_pch_all <- color_pch[match(as.character(in_phy$reference), names(color_pch))]
  }
  if (!is.null(transparent_pch)) {
    color_pch_all <- gplots::col2hex(color_pch_all)
    color_pch_all <- paste0(color_pch_all, transparent_pch)
  }
  for (i in unique(in_phy$mrca_node_number)) {
    rowsies <- in_phy$mrca_node_number == i
    x_min <- min(x_ages[rowsies])
    x_max <- max(x_ages[rowsies])
    # we just need one element for y position of nodes:
    y_pos <- y_nodes[rowsies][1]
    graphics::segments(x0 = x_min,
             y0 = y_pos,
             x1 = x_max,
             y1 = y_pos,
             col = color_bars,
             lty = lty_bars,
             lwd = lwd_bars)
  }
  graphics::points(x_ages,
         y_nodes,
         col = color_pch_all,
         pch = pch,
         cex = cex_pch)
  # add a time axis
  ## choose an axis plot_type
  match.arg(arg = plot_type, choices = c("phyloch", "ape"))
  ## get the geologic timescale
  if (is.null(geologic_timescale) | "strat2012" %in% geologic_timescale) {
    utils::data("strat2012", package = "phyloch")
    geologic_timescale <- strat2012
  }
  if ("ape" %in% plot_type) {
    # TODO: fix axis not plotting all the way through the root
    # See issue https://github.com/phylotastic/datelifeplot/issues/1
    ape::axisPhylo(side = 1, line = pos_axis, cex.axis = cex_axis)
    axisChrono(side = 1, unit = NULL, line = pos_axis, cex.axis = cex_axis)
  } else { # if ("phyloch" %in% plot_type) {
    axisGeo(GTS = geologic_timescale,
            unit = geologic_unit,
            col = c("gray80", "white"),
            gridcol = c("gray80", "white"),
            cex = cex_axis,
            gridty = "twodash")
  }
  # add a label to the axis
  if (!is.null(units)) {
    graphics::mtext(axis_label,
                    cex = cex_axislabel,
                    side = 1,
                    font = 2,
                    line = (omi1-0.2)/0.2,
                    outer = FALSE,
                    at = max(lastPP$xx) * center_axislabel)# centering of the time axis label
  }
  # add a title to the plot
  if (!is.null(title)) {
    titlei <- wrap_string_to_plot(string = title, max_cex = cex_title, whole = FALSE)
    graphics::mtext(text = titlei$wrapped, outer = TRUE,
      cex = titlei$string_cex, font = titlei$string_font, line = pos_title)
  }
  if (add_legend) {
    x_legend <- ifelse(is.null(x_legend), -time_depth*0.5, x_legend)
    y_legend <- ifelse(is.null(y_legend), max(y_nodes), y_legend)

    if (x_legend <= 0) {
      graphics::par(xpd=TRUE) # so it's clipped in the outer margin
    }
    if (y_legend > max(y_nodes)) {
      graphics::par(xpd=NA) # so it's clipped in the upper margin
    }
    message("Current legend x co-ordinate is set to ", x_legend)
    message("And y co-ordinate is set to ", y_legend)
    graphics::legend(x = x_legend,
           y = y_legend,
           legend = names(color_pch),
           pch = 19,
           col = color_pch,
           cex = cex_legend)
  }
}
