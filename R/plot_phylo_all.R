#' Plot all Chronograms with Study Titles and Geochronological Axis
#'
#' @param chronograms A list of chronograms as `multiPhylo` or as a plain `list` object.
#' @inheritParams subset_trees
#' @inheritParams plot_phylo
#' @param individually Boolean indicating if chronograms should be plotted one by one
#' or appended to the same file.
#' @inheritDotParams ape::plot.phylo
#' @param folder_name A character string indicating the name of the folder to write
#' the chronogram plot files to. Default to "phylo_all". Only relevant if `write = "png"` or `write = "pdf"`.
#' @param file_name A character string indicating a file name to write the chronogram
#' plots to. Relevant if `write = "png"` or `write = "pdf"`.
#' @param max_depth A numeric vector of length 1, indicating the upper limit of
#'   the time scale on the x axis to be used on all plots. If none is provided, it is estimated
#'   by getting the time depth of the oldest chronogram, adding 5 and rounding to the
#'   closest 10. See details for more.
#' @details
#' Currently, max_depth is obtained by default with `round(max(sapply(chronograms, function(x) max(ape::branching.times(x)))) + 5, digits = -1)`.
#' @export
plot_phylo_all <- function(chronograms,
                           include = TRUE,
                           individually = TRUE,
                           folder_name = "phylo_all",
                           max_depth, # time_depth in plot_phylo
                           plot_type = "phyloch",
                           time_axis = TRUE,
                           mai1, mai2, mai3, mai4,
                           omi1, omi2, omi3, omi4,
                           plot_height, plot_width,
                           write = "no",
                           file_name = "chronogram",
                           geologic_timescale = "strat2012",
                           geologic_unit = "period",
                           axis_label = "Time (MYA)",
                           cex_tiplabels = graphics::par("cex"),
                           cex_axislabel = graphics::par("cex"),
                           cex_axis = graphics::par("cex"),
                           cex_title = graphics::par("cex"),
                           pos_title = 1,
                           pos_axis = 1,
                           center_axislabel = 0.5,
                           ...) {
  chronograms <- subset_trees(chronograms, include = include)
  # in case there is just one tree in chronograms
  if (any("tip.label" %in% names(chronograms))) {
    chronograms <- list(chronograms)
  }
  if (missing(max_depth)) {
    max_depth <- round(max(sapply(chronograms,
                                  function(x) max(ape::branching.times(x)))) + 5, digits = -1)
    # in case we used it for a datelife object:
    # if(isTRUE(all.equal(round(sapply(chronograms,
    # function(x) max(ape::branching.times(x))), digits = 3))))
    # max_depth <- round(max(summarize_datelife_result(datelife_result = get.filtered.results(),
    #             summary_format = "mrca")) + 5)

  }
  if (!is.numeric(max_depth)) {
    stop("'max_depth' argument must be 'numeric'. It currently is ", max_depth)
  }
  # to plot all chronograms from the same depth, add a root:
  # using lapply was working, until it stopped working, so I replaced with a for loop
  # chronograms <- lapply(chronograms, function(x) {
  #   x$root.edge <- max_depth - max(ape::branching.times(x))
  # })
  for (i in seq(chronograms)) {
    chronograms[[i]]$root.edge <- max_depth - max(ape::branching.times(chronograms[[i]]))
  }
  ##############################################################################
  ##############################################################################
  # Getting mai4 data
  ##############################################################################
  ##############################################################################
  tip_labels <- unique(unlist(sapply(chronograms, "[", "tip.label")))
  ind <- which.max(nchar(tip_labels))
  mai4_in <- graphics::strwidth(s = tip_labels[ind],
                             units = "inches",
                             cex = cex_tiplabels,
                             font = 3)
  if(missing(mai4)) {
    mai4 <- mai4_in
  }
  # if(any(lapply(chronograms, ape::Ntip) > 3))
  # png("~/tmp/axisgeo.png", units = "in")
  if (!(write %in% c("png", "pdf"))) {
  # if (!grDevices::devAskNewPage()
  # && !names(grDevices::dev.cur()) %in% c("pdf", "postscript")) {
      grDevices::devAskNewPage(TRUE)
      # dev.size("px")
      on.exit(grDevices::devAskNewPage(FALSE))
  } else {
      if (!dir.exists(folder_name)) {
        dir.create(path = gsub("\\.png$|\\.pdf$", "", folder_name))
      }
  }

  if (length(file_name) != length(chronograms)) {
    file_name <- paste0(folder_name,
                        "/",
                        gsub("\\.png$|\\.pdf$", "", file_name),
                        "_")
    file_name <- paste0(file_name,
                        seq(chronograms),
                        ".",
                        write)
  } else {
    file_name <- paste0(folder_name,
                        "/",
                        gsub("\\.png$|\\.pdf$", "", file_name),
                        ".",
                        write)
  }
  if (length(plot_height) != length(chronograms)) {
    plot_height <- rep(plot_height, length(chronograms))
  }
  if (length(plot_width) != length(chronograms)) {
    plot_width <- rep(plot_width, length(chronograms))
  }
  for (i in seq(chronograms)) {
    plot_phylo(chronogram = chronograms[[i]],
               title = names(chronograms)[i],
               time_depth = max_depth,
               plot_type = plot_type,
               time_axis = TRUE,
               mai1 = mai1,
               mai2 = mai2,
               mai3 = mai3,
               mai4 = mai4,
               omi1 = omi1,
               omi2 = omi2,
               omi3 = omi3,
               omi4 = omi4,
               plot_height = plot_height[i],
               plot_width = plot_width[i],
               write = write,
               file_name = file_name[i],
               geologic_timescale = NULL,
               geologic_unit = geologic_unit,
               axis_label = axis_label,
               cex_tiplabels = cex_tiplabels,
               cex_axislabel = cex_axislabel,
               cex_axis = cex_axis,
               cex_title = cex_title,
               pos_title = pos_title,
               pos_axis = pos_axis,
               center_axislabel = center_axislabel,
               ...)
  }
  # getting an "unrecoverable" error with merge PDF:
  # if (!individually) {
  #   # install_github("trinker/plotflow")
  # if we decide to use this, we should add plotflow functions in datelife package so we don't have to add it to description...
  #   plotflow:::mergePDF(
  #       in.file= paste(file.path(gsub("\\.png$|\\.pdf$", "", file), dir(gsub("\\.png$|\\.pdf$", "", file)))),
  #       file= paste0(gsub("\\.png$|\\.pdf$", "", file), ".", write)
  #       # file= "merged.pdf"
  #   )
  # }
}
