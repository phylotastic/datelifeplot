# @details The node position on the y axis is given by pp$yy
# pp$yy also stores the positions of the tips
# the y coordinate of nodes start at pp$yy[ape::Ntip + 1]

#' Plot a Chronogram and Add Age Data to Corresponding Nodes
#' @description `plot_node_ages` plots a given chronogram and adds age
#'   data [points()] given in `node_ages` and `summary_ages` to the corresponding nodes.
#'
#' @param chronogram A `phylo` object with branch length proportional to time.
#' @param matched_ages An output of [datelife:::summary.matchedCalibrations()] with age data matched to nodes in `chronogram`.
#' @param pch_color A named vector of colors. Names must correspond to names
#'  in `matched_ages$in_phy$references`.
#'  If vector is not named and length > 1, colors will be recycled.
#' @param pch_transparency A numeric value ranging from 10-99, indicating
#' color transparency for color defined by `pch_color`. Default to `NULL`, no transparency.
#' @param pch A numeric vector indicating the symbol for calibration age points.
#'   See [graphics::par()] for options. Defaults to 20 = "bullet circle".
#' @param pch_cex A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of symbols defined by `pch`. Default to 1.
#' @param bars_color A character vector of one element indicating the color for
#'  node age distribution bars. Defaults to "#80808050", which is hex for gray
#' ("#808080") with 50\% transparency.
#' @param bars_lty A numeric vector or character string indicating the plotting
#'   line type for bars. See [graphics::par()] for options. Default to "solid".
#' @param bars_lwd A numeric vector indicating the line width for age distribution bars.
#'   See [graphics::par()] for options. Default to 7.
#' @param add_legend Default to `TRUE`, adds a legend to the left of the plot.
#' @param legend_cex A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of legend. Default to one half the size of the axis label, `cex_axislabel * 0.5`.
#' @param legend_x the x co-ordinate to be used to position the legend on the left side of the plot.
#' @param legend_y the y co-ordinate to be used to position the legend on the left side of the plot.
#' @param legend_box Default to `TRUE`, adds a box around the legend.
#' @inheritParams plot_phylo
#' @inheritDotParams ape::plot.phylo
#' @importFrom ape .PlotPhyloEnv
#'@details Plot are margin sizes as defined by [graphics::par()$mai] and [graphics::par()$omi]
#' are overruled within the function. To modify them you have to use the arguments
#' `mai1`, `mai2`, `mai3` and `mai4`, and omi1, omi2, omi3 and omi4.
#' @export
plot_node_ages2 <- function(chronogram,
                           matched_ages,
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
                           pch = 20,
                           pch_color,
                           pch_transparency = NULL,
                           pch_cex = graphics::par("cex"),
                           bars_color = "#80808050",
                           bars_lty = "solid",
                           bars_lwd = 7,
                           add_legend = TRUE,
                           legend_cex = cex_axislabel*0.5,
                           legend_x = NULL,
                           legend_y = NULL,
                           legend_box = TRUE,
                           ...) {
  #
  if (missing(matched_ages)) {
    stop("argument 'matched_ages' is missing with no default.")
  }
  # get calibrations that are in_phy only
  # creates a list of data.frames:
  in_phy <- lapply(matched_ages, "[[", "in_phy")
  # obtain max x lim from calibration ages, chronogram and chronogram + root:
  phylo_length <- max(ape::branching.times(chronogram))
  max_ages <- unlist(lapply(in_phy, function(x) x[,"MaxAge"]))
  min_ages <- unlist(lapply(in_phy, function(x) x[,"MinAge"]))
  max_calibration  <- max(c(max_ages, min_ages))
  if (is.null(time_depth)) {
    if (is.null(chronogram$root.edge)) {
      time_depth <- round(max(phylo_length, max_calibration) * 1.2, digits = -1)
    } else {
      time_depth <- round(max(phylo_length + chronogram$root.edge, max_calibration), digits = -1)
    }
  }
  ############################################################################
  # start chunk that can be replaced by plot_chronogram
  # plot_chronogram(chronogram,
  #                       title = "Chronogram",
  #                       time_depth = NULL,
  #                       plot_type = "phyloch",
  #                       time_axis = TRUE,
  #                       mai1, mai2, mai3, mai4,
  #                       omi1, omi2, omi3, omi4,
  #                       plot_height, plot_width,
  #                       geologic_timescale = "strat2012",
  #                       geologic_unit = "period",
  #                       cex_tiplabels = graphics::par("cex"), # inherits param from plot_node_ages
  #                       cex_axislabel = graphics::par("cex"),
  #                       cex_axis = graphics::par("cex"),
  #                       cex_title = graphics::par("cex"),
  #                       pos_title = 1,
  #                       pos_axis = 1,
  #                       center_axislabel = 0.5,
  #                       axis_label = "Time (MYA)")
  ############################################################################
  ############################################################################
  # add a root to the chronogram (or extend the root) to plot from the specified time_depth
  chronogram$root.edge <- time_depth - max(ape::branching.times(chronogram))
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
    ind <- which.max(nchar(chronogram$tip.label))
    mai4 <- graphics::strwidth(s = chronogram$tip.label[ind],
                               units = "inches",
                               cex = cex_tiplabels,
                               font = 3)
  }
  if (missing(omi1)) {
    pho <- phylo_height_omi(phy = chronogram)
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
                mai = c(mai1, mai2, mai3, mai4),
                omi = c(omi1, omi2, omi3, omi4))

  # plot.phylo FALSE
  ape::plot.phylo(chronogram,
                  cex = cex_tiplabels,
                  label.offset = 0.5,
                  x.lim = c(0, time_depth),
                  root.edge = TRUE,
                  plot = FALSE,
                  ...)
  # get plotting x.lim
  firstPP_x.lim <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2]
  ape::plot.phylo(chronogram,
                  cex = cex_tiplabels, #edge.width = 2,
                  label.offset = 0.5,
                  x.lim = c(0, firstPP_x.lim),
                  root.edge = TRUE,
                  plot = TRUE, ...)

  # get recorded plotting parameters
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  # END chunk that can be replaced by plot_chronogram
  ############################################################################

  ############################################################################
  ############################################################################
  # plot gray bars
  ############################################################################
  ############################################################################
  # dplyr::bind_reows only works if class data.frame
  for (i in seq(in_phy)) {
    class(in_phy[[i]]) <- "data.frame"
  }
  combined <- dplyr::bind_rows(in_phy)
  # we use lastPP$xx positions to get plot x position for bar age range
  ages_equal <- combined$MinAge == combined$MaxAge
  if (all(ages_equal)) {
    # We need the one max x position in lastPP$xx with max(lastPP$xx)
    # Then we substract the nide ages to get actual x positions for
    # all node ages on the plot:
    x_ages <- max(lastPP$xx) - combined$MinAge
  }
  # Next, we use lastPP$yy positions to get plotting y position of nodes
  if (all(combined$mrca_node_number > lastPP$Ntip)) {
    # case in which mrca_node_numbers start at Ntip + 1
    # we can directly use mrca_node_numbers as index for lastPP$yy
    # it generates a vector of the correct length
    y_ages <- lastPP$yy[combined$mrca_node_number]
  }
  for (i in unique(combined$mrca_node_number)) {
    rowsies <- combined$mrca_node_number == i
    x_min <- min(x_ages[rowsies])
    x_max <- max(x_ages[rowsies])
    # we just need one element for y position of nodes:
    y_pos <- y_ages[rowsies][1]
    graphics::segments(x0 = x_min,
             y0 = y_pos,
             x1 = x_max,
             y1 = y_pos,
             col = bars_color,
             lty = bars_lty,
             lwd = bars_lwd)
  }
  ############################################################################
  ############################################################################
  # plot points
  ############################################################################
  ############################################################################
  legend_color <- legend_pch <- c()
  for (matched in names(in_phy)) {
    # we use lastPP$xx positions to get plot x position for node ages
    ages_equal <- in_phy[[matched]]$MinAge == in_phy[[matched]]$MaxAge
    if (all(ages_equal)) {
      # We need the one max x position in lastPP$xx with max(lastPP$xx)
      # Then we substract the node ages to get actual x positions for
      # all node ages on the plot:
      x_ages <- max(lastPP$xx) - in_phy[[matched]]$MinAge
    }
    # Next, we use lastPP$yy positions to get plotting y position of nodes
    if (all(in_phy[[matched]]$mrca_node_number > lastPP$Ntip)) {
      # case in which mrca_node_numbers start at Ntip + 1
      # we can directly use mrca_node_numbers as index for lastPP$yy
      # it generates a vector of the correct length
      y_ages <- lastPP$yy[in_phy[[matched]]$mrca_node_number]
    }
    # Use study references to color points
    if (missing(pch_color)) {
      color_pch_all <- in_phy[[matched]]$reference
      print(names(color_pch_all ))
    } else {
      mm <- match(as.character(in_phy[[matched]]$reference), names(pch_color[[matched]]))
      color_pch_all <- pch_color[[matched]][mm]
      print(names(color_pch_all ))
    }
    # assign transparency
    if (!is.null(pch_transparency)) {
      color_pch_all <- gplots::col2hex(color_pch_all)
      color_pch_all <- paste0(color_pch_all, pch_transparency)
    }
    legend_color <- c(legend_color, color_pch_all)
    # plot gray bars:
    for (i in unique(in_phy[[matched]]$mrca_node_number)) {
      rowsies <- in_phy[[matched]]$mrca_node_number == i
      x_min <- min(x_ages[rowsies])
      x_max <- max(x_ages[rowsies])
      # we just need one element for y position of nodes:
      y_pos <- y_ages[rowsies][1]
      graphics::segments(x0 = x_min,
               y0 = y_pos,
               x1 = x_max,
               y1 = y_pos,
               col = bars_color,
               lty = bars_lty,
               lwd = bars_lwd)
    }
    # plot in_phy ages
    graphics::points(x_ages,
           y_ages,
           col = color_pch_all,
           pch = pch[[matched]],
           cex = pch_cex)
   xx <- rep(pch[[matched]], length(pch_color[[matched]]))
   legend_pch <- c(legend_pch, xx)
}
  if (add_legend) {
    legend_text <- unlist(sapply(pch_color, names))
    if (length(legend_text) == 0) {
      stop("argument 'pch_color' must have names corresponding to legend.")
    }
  }
  ############################################################################
  ############################################################################
  # add a time axis
  ############################################################################
  ############################################################################
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
  if (!is.null(axis_label)) {
    graphics::mtext(axis_label,
                    cex = cex_axislabel,
                    side = 1,
                    font = 2,
                    line = (omi1-0.2)/0.2,
                    outer = FALSE,
                    at = max(lastPP$xx) * center_axislabel)# centering of the time axis label
  }
  ############################################################################
  ############################################################################
  # add a title to the plot
  ############################################################################
  ############################################################################
  if (!is.null(title)) {
    titlew <- wrap_string_to_plot(string = title, max_cex = cex_title, whole = FALSE)
    graphics::mtext(text = titlew$wrapped, outer = TRUE,
      cex = titlew$string_cex, font = titlew$string_font, line = pos_title)
  }
  ############################################################################
  ############################################################################
  # add a legend
  ############################################################################
  ############################################################################
  if (add_legend) {
    legend_x <- ifelse(is.null(legend_x), -time_depth*0.5, legend_x)
    legend_y <- ifelse(is.null(legend_y), max(y_ages), legend_y)

    if (legend_x <= 0) {
      graphics::par(xpd=TRUE) # so it's clipped in the outer margin
    }
    if (legend_y > max(y_ages)) {
      graphics::par(xpd=NA) # so it's clipped in the upper margin
    }
    message("Current legend x co-ordinate is set to ", legend_x)
    message("And y co-ordinate is set to ", legend_y)
    graphics::legend(x = legend_x,
           y = legend_y,
           legend = unlist(legend_text),
           pch = legend_pch,
           col = unlist(legend_color),
           cex = legend_cex,
           bty = ifelse(legend_box, "o", "n")
         )
  }
}
