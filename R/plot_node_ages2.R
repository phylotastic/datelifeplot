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
#' @param legend_pch A numeric vector indicating the symbol to use in legend.
#'   See [graphics::par()] for options. Defaults to 20 = "bullet circle".
#' @param legend_text A character vector indicating the text to use as legend.
#' @param legend_color A character vector indicating the colors to use in legend.
#' @param legend_title A character vector indicating the title to use in legend.
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
                           legend_x,
                           legend_y,
                           legend_box = TRUE,
                           legend_pch,
                           legend_text,
                           legend_color,
                           legend_title,
                           ...) {
  #
  if (missing(matched_ages)) {
    stop("argument 'matched_ages' is missing with no default.")
  }
  # get calibrations that are in_phy only
  # create a list of data.frames:
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
  # sometimes, nodes have only one age point, or ages that do not go over the node
  # to have gray bars connected to the nodes, we have to include their age in the range
  # first, get chronogram node ages :
  branching_times <- ape::branching.times(chronogram)
  if (is.null(names(branching_times))) {
    stop("nodes in chronogram need to be named.")
  } else {
    message("'chronogram' node names used to extract node numbers are \n",
            paste(head(names(branching_times)), collapse = ", "),
            " ... ", paste(tail(names(branching_times)), collapse = ", "))
  }
  chronogram_node_numbers <- as.numeric(sub("n", "", names(branching_times)))
  # create a chronogram node ages data frame to combine:
  times_data_frame <- data.frame(MRCA = names(branching_times),
                                 MaxAge = as.numeric(branching_times),
                                 MinAge = as.numeric(branching_times),
                                 taxonA = rep("taxonA", length(branching_times)),
                                 taxonB = rep("taxonB", length(branching_times)),
                                 reference = rep("chronogram", length(branching_times)),
                                 mrca_node_number = chronogram_node_numbers,
                                 mrca_node_name = names(branching_times),
                                 nodeAge = as.numeric(branching_times))
  # Combine calibrations and summary ages:
  # dplyr::bind_rows only works for class data.frame, so convert ages to data.frame:
  for (i in seq(in_phy)) {
    class(in_phy[[i]]) <- "data.frame"
  }
  ## combine the data.frames into a single one:
  combined <- dplyr::bind_rows(c(in_phy, list(times_data_frame)))
  # homogenize node numbers to be Ntip + 1
  below_ntip <- combined$mrca_node_number <= lastPP$Ntip
  fixed <- combined$mrca_node_number[below_ntip] +  lastPP$Ntip
  combined$mrca_node_number[below_ntip] <- fixed
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
  } else {
    stop("There is somethong wrong with mrca node numbers in 'matched_ages'")
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
  legend_color_in <- legend_pch_in <- pch_color_in <- vector(mode = "list")
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
      # choose colors at random
      color_pch_all <- in_phy[[matched]]$reference
    } else {
      # choose colors from pch_color, matching names:
      mm <- match(as.character(in_phy[[matched]]$reference), names(pch_color[[matched]]))
      color_pch_all <- pch_color[[matched]][mm]
    }
    # assign transparency
    if (!is.null(pch_transparency)) {
      color_pch_all <- gplots::col2hex(color_pch_all)
      color_pch_all <- paste0(color_pch_all, pch_transparency)
    }
    pch_color_in <- c(pch_color_in, list(color_pch_all))
    legend_color_in <- c(legend_color_in, list(color_pch_all))
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
    # plot points of in_phy ages
    graphics::points(x_ages,
           y_ages,
           col = color_pch_all,
           pch = pch[[matched]],
           cex = pch_cex)
   # get pch symbols for legend:
   xx <- rep(pch[[matched]], length(pch_color[[matched]]))
   legend_pch_in <- c(legend_pch_in, list(xx))
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
                    at = max(lastPP$xx) * center_axislabel)  # centering of the time axis label
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
  if (any(add_legend)) {
    if (missing(legend_title)) {
      legend_title <- paste("Data set", seq(matched_ages))
    }
    if (missing(legend_x)) {
      legend_x <- -time_depth * 0.5
    }
    if (missing(legend_y)) {
      legend_y <- max(y_ages)
    }
    if (length(legend_x) != length(matched_ages)) {
      legend_x <- rep(legend_x, length(matched_ages))
    }
    if (length(legend_y) != length(matched_ages)) {
      legend_y <- rep(legend_y, length(matched_ages))
    }
    if (missing(legend_cex)) {
      legend_cex <- cex_axislabel * 0.5
    }
    if (length(legend_cex) != length(matched_ages)) {
      legend_cex <- rep(legend_cex, length(matched_ages))
    }
    if (length(legend_box) != length(matched_ages)) {
      legend_box <- rep(legend_box, length(matched_ages))
    }
    # add a legend for each data set in matched_ages
    for (i in seq(matched_ages)) {
      # determine text for legend:
      if (missing(legend_text)) {
        legend_text_i <- names(pch_color_in[[i]])
        if (length(legend_text_i) == 0) {
          message("There is no text for legend.\n",
          "You can provide one in the 'legend_text' argument or as names in 'pch_color' argument.")
          warning("Legend was not added.")
          break()
        }
      } else {
        legend_text_i <- legend_text[[i]]
      }
      if (inherits(legend_text_i, "list")) {
        legend_text_i <- unlist(legend_text_i)
      }
      # determine legend cex:
      legend_cex_i <- legend_cex[i]
      # determine legend pch:
      legend_pch_i <- unlist(ifelse(missing(legend_pch),
                             legend_pch_in[i],
                             legend_pch[i]))
      # determine legend color:
      legend_color_i <- unlist(ifelse(missing(legend_color),
                               legend_color_in[i],
                               legend_color[i]))
      # determine x and y position of legend:
      legend_x_i <- legend_x[i]
      legend_y_i <- legend_y[i]
      # legend box
      legend_box_i <- ifelse(legend_box[i], "o", "n")
      # legend title
      legend_title_i <- legend_title[i]
      # using NULL for "yes" argument of ifelse() errored as follows:
      # Error in ans[ypos] <- rep(yes, length.out = len)[ypos] :
      #   replacement has length zero
      # In addition: Warning message:
      # In rep(yes, length.out = len) : 'x' is NULL so the result will be NULL
      # determine xpd par() so legend is not hidden. We have two cases:
      if (legend_x_i <= 0) { ## if legend goes in the left margin
        graphics::par(xpd = TRUE) # xpd = TRUE clips the left margin
      }
      if (legend_y_i > max(y_ages)) { ## if legend goes in the upper margin
        graphics::par(xpd = NA) # xpd = NA clips the upper margin
      }
      message("legend ", i, " x co-ordinate is set to ", legend_x_i)
      message("And ", i, " y co-ordinate is set to ", legend_y_i)
      # plot actual legend:
      graphics::legend(x = legend_x_i,
             y = legend_y_i,
             legend = legend_text_i,
             pch = legend_pch_i,
             col = legend_color_i,
             cex = legend_cex_i,
             bty = legend_box_i,
             title = legend_title_i)
    }
  }
}
