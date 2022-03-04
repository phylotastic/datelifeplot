#' Plot a single chronogram with a title and a geochronological axis
#'
#' @details [plot_phylo()] uses different plotting functions to generate a plot
#' of a given chronogram with a time axis representing geologic time.
#'
#' @param chronogram A chronogram either as a newick character string or as a
#' `phylo` object with branch length proportional to time.
#' @param title A character string providing the title for the plot.
#' @param time_depth A numeric vector indicating the upper limit for the time scale on the x axis.
#' @param plot_type A character vector of length one indicating the type of
#' plot to generate. Options are:
#' \describe{
#' 	\item{"ape"}{It uses the functions [ape::plot.phylo()] and [ape::axisPhylo()].}
#' 	\item{"phyloch"}{It uses the functions [ape::plot.phylo()] and [phyloch::axisGeo()].}
#' 	\item{"strap"}{It uses the function [strap::geoscalePhylo()] from the package [strap].}
#' 	\item{"phytools"}{Not implemented yet. It will use functions from the package [phytools]}
#' 	}
#' @param time_axis Default to `TRUE`, add a time axis to the tree plot.
#' @param mai1,mai2,mai3,mai4 A numeric value indicating internal plot margin sizes
#' in inches. `mai4` is particularly important as it indicates the space needed for
#' plotting whole tip labels on the right margin of the plot.
#' @param omi1,omi2,omi3,omi4 A numeric value indicating outter plot margin sizes in inches.
#' @param plot_height,plot_width A numeric value indicating height and width for the plot.
#' @param geologic_timescale A dataframe of geochronological limits.
#' @param geologic_unit A character vector used to select geological time units that shall be displayed. When using \code{gradstein04}, \code{"eon"}, \code{"era"}, \code{"period"}, \code{"epoch"}, and \code{"stage"} are available.
#' @param cex_axis A numeric indicating character expansion for the axis. Default to value given by `graphics::par("cex")`.
#' @param cex_title A numeric indicating character expansion for the title. Default to value given by `graphics::par("cex")`.
#' @param cex_axislabel A numeric indicating character expansion for the time axis label. Default to value given by `graphics::par("cex")`.
#' @param pos_title Indicates the line position of the title. Default to 1.
#' @param pos_axis Indicates the line position of the axis. Default to 1.
#' @param center_axislabel A numeric indicating centering position of time axis label. Default to 0.5.
#' @param axis_label A character string used to provide information on time units. Defaults
#' to "Time (MYA)" (time in million years ago). If NULL, time units are not added.
#' @param cex_tiplabels A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of tip labels in `chronogram`. Default to value given by `graphics::par("cex")`.
#' @inheritParams phyloch::axisGeo
#' @inheritDotParams ape::plot.phylo
#' @importFrom ape .PlotPhyloEnv
#' @export
# enhance: examples of axis_types!
plot_chronogram <- function(chronogram,
                       title = "Chronogram",
                       time_depth = NULL,
                       plot_type = "phyloch",
                       time_axis = TRUE,
                       mai1, mai2, mai3, mai4,
                       omi1, omi2, omi3, omi4,
                       plot_height, plot_width,
                       geologic_timescale = "strat2012",
                       geologic_unit = "period",
                       cex_tiplabels = graphics::par("cex"), # inherits param from plot_node_ages
                       cex_axislabel = graphics::par("cex"),
                       cex_axis = graphics::par("cex"),
                       cex_title = graphics::par("cex"),
                       pos_title = 1,
                       pos_axis = 1,
                       center_axislabel = 0.5,
                       axis_label = "Time (MYA)",
                       #pos_axislabel = ,
                       ...){
  #
  print(time_depth)
  if (!inherits(chronogram, "phylo")) {
      if (inherits(chronogram, "multiPhylo")) {
        message("'chronogram' is a 'multiPhylo' object. Only first chronogram will be plotted.")
        chronogram <- chronogram[[1]]
      } else {
        message("'chronogram' must be a 'phylo' or 'multiPhylo' object.")
        return(NA)
      }
  }
  if (is.null(chronogram$edge.length)) {
    message("'chronogram' has no edge lengths; a geologic time axis can not be plotted.")
    return(NA)
  }
  phylo_length <- max(ape::branching.times(chronogram))
  if (is.null(time_depth)) {
    if (is.null(chronogram$root.edge)) {
      time_depth <- round(phylo_length*1.2, digits = -1)
    } else {
      time_depth <- round(phylo_length + chronogram$root.edge, digits = -1)
    }
  }
  # add a root to the chronogram (or extend the root) to plot from the specified time_depth
  chronogram$root.edge <- time_depth - max(ape::branching.times(phy))
  match.arg(arg = plot_type, choices = c("phyloch", "strap", "phytools", "ape"))
  if (is.null(geologic_timescale) | "strat2012" %in% geologic_timescale) {
    utils::data("strat2012", package = "phyloch")
    geologic_timescale <- strat2012
  }
  if (missing(mai1)) {
    mai1 <- 0
  }
  if (missing(mai2)) {
    mai2 <- 0
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
  pho <- phylo_height_omi(phy = chronogram)
  message("Recommended plot area height is ", pho$height)
  if (missing(plot_height)) {
    plot_height <- pho$height
  }
  if (missing(omi1)) {
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
  # plot_chronogram.phylo(chronograms[[i]], cex = 1.5, edge.width = 2, label.offset = 0.5,
    # x.lim = c(0, max_depth), root.edge = TRUE, root.edge.color = "white")
  # graphics::par(xpd = FALSE)
  print(paste0("chronogram$root.edge = ", chronogram$root.edge))
  if (plot_type %in% c("ape", "phyloch")) {
    ape::plot.phylo(chronogram,
                    cex = cex_tiplabels, #edge.width = 2,
                    label.offset = 0.5,
                    x.lim = c(0, time_depth),
                    root.edge = TRUE,
                    plot = FALSE, ...)  #
    firstPP_x.lim <- get("last_plot.phylo", envir = .PlotPhyloEnv)$x.lim[2]
    ape::plot.phylo(chronogram,
                    cex = cex_tiplabels, #edge.width = 2,
                    label.offset = 0.5,
                    x.lim = c(0, firstPP_x.lim),
                    root.edge = TRUE,
                    plot = TRUE, ...)  #
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if(time_axis) {
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
    }
    # center_axislabel <- 0.5
  }
  if ("phytools" %in% plot_type) {
    # TODO
    message("Plotting a geologic time axis with phytools is not supported yet.")
  }
  if ("strap" %in% plot_type) {
    chronogram$root.time <- phylo_length
    strap::geoscalePhylo(tree = chronogram,
                         x.lim = c(0, phylo_length),
                         cex.tip = 0.7,
                         cex.ts = 0.7,
                         cex.age = 0.7,
                         width = 4,
                         tick.scale = 15,
    # create boxes using the smallest geologic unit in argument "geologic_unit":
                         boxes = geologic_unit[length(geologic_unit)],
                         erotate = 90,
                         quat.rm = TRUE,
                         units = axis_label)
    # center_axislabel <- 1
  }
  # add a label to the axis
  if (time_axis & !is.null(axis_label)) {
    graphics::mtext(axis_label,
                    cex = cex_axislabel,
                    side = 1,
                    font = 2,
                    line = (omi1-0.2)/0.2,
                    outer = FALSE,
                    at = max(lastPP$xx) * center_axislabel # centering of the time axis label
    )
  }
  # add a title to the plot
  if (!is.null(title)) {
    titlei <- wrap_string_to_plot(string = title, max_cex = cex_title, whole = FALSE)
    graphics::mtext(text = titlei$wrapped, outer = TRUE,
      cex = titlei$string_cex, font = titlei$string_font, line = pos_title)
  }
}
