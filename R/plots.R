#' A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.
#'
#' @name strat2012
#' @docType data
#' @format A data frame
#' \describe{
#'   \item{eon}{Eon name.}
#'   \item{era}{Era name.}
#'   \item{period}{Period name.}
#'   \item{epoch}{Epoch name.}
#'   \item{stage}{Stage name.}
#'   \item{MA}{Estimated boundary age for the associated interval.}
#'   \item{error}{Estimated errors associated with each age estimate.}
#'   \item{GSSP}{Binary response denoting whether the age estimate is defined by a basal Global Standard Section and Point}
#' }
#' @source \url{https://github.com/fmichonneau/phyloch}
#' @keywords plot gradstein geochronology
#' @details
#' See phyloch package data description.
#' Enhance: There are a couple errors in the original strat2012 object that are not too meaningful but could be improved.
#' Generated with utils::data("strat2012", package = "phyloch")
#' use_data(strat2012)
"strat2012"

#' subset chronograms for plotting densitree plots and phylo_all plots
#' @param chronograms A list of chronograms as multiPhylo or as a plain list object.
#' @param include Boolean or numeric vector. Default to TRUE, keep all chronograms
#' in chronograms. If FALSE, exclude chronograms with only two tips. If numeric, it is used
#' as indices to subset chronograms object.
subset_trees <- function(chronograms, include = TRUE){
  if(is.numeric(include)){  #if it is numeric
    include <- round(include)
    include <- include[which(include <= length(chronograms))]
    include <- sort(unique(include))
    chronograms <- chronograms[include]
  }
  if(is.logical(include) & !is.na(include)){
      if (!include){
      chronograms <- chronograms[sapply(chronograms, function (x) ape::Ntip(x)) > 2]
    }
  }
  chronograms
}

#' Get a densiTree plot from a set of opentree_chronograms
#' @details If [phangorn::densiTree()] throws an error, [plot_densitree()] chooses
#' the chronogram with the most tips as consensus (using [datelife::get_biggest_phylo()])
#' we found that [phangorn::densiTree()] errors commonly from failing to do a consensus tree.
#' @param chronograms A list of chronograms as multiPhylo or as a plain list object.
#' @inheritParams subset_trees
#' @inheritDotParams phangorn::densiTree -x -consensus
#' @export
plot_densitree <- function(chronograms, include = TRUE, ...) {
  chronograms <- subset_trees(chronograms, include = include)
  # for densitree plot it does not really matter if chronograms have the same length
  # max_depth <- round(max(sapply(chronograms, function(x) max(ape::branching.times(x)))) + 5, digits = -1)
  # # max_depth <- round(max(summarize_datelife_result(datelife_result = get.filtered.results(),
  # #             summary_format = "mrca")) + 5)  # in case we used it for a datelife object
  #
  # # plot all chronograms from the same depth:
  # chronograms <- lapply(chronograms, function(x) {
  #  x$root.edge <- max_depth - max(ape::branching.times(x))
  #  x
  # })
  class(chronograms) <- "multiPhylo"
  # if we use biggest phylo as consensus for all, some data set are not plotted correctly
  # biggest_phylo <- datelife::get_biggest_multiphylo(chronograms = chronograms)
  # try(phangorn::densiTree(x = chronograms, consensus = biggest_phylo, ...))
  tryCatch(phangorn::densiTree(x = chronograms, ...),
  error = function(e) {
    biggest_phylo <- datelife::get_biggest_multiphylo(trees = chronograms)
    try(phangorn::densiTree(x = chronograms, consensus = biggest_phylo, ...))
  })
}

#' get the outer margin of a graphics device from a number of tips
#' @param phy `phylo` object to be plotted
phylo_height_omi <- function(phy){
  tipnum <- ape::Ntip(phy)
  if(tipnum > 10){
    hei <- 50 + (30 * tipnum)
  } else {
    hei <- 300
  }
  if(tipnum == 2){
    omi1 <- 2
  } else if(tipnum == 3){
    omi1 <- 1.5
  } else if(tipnum == 4){
    omi1 <- 1.2
  } else if (tipnum >= 5 & tipnum <= 7){
    # omi1 <- 3
    omi1 <- 1 # this works good for pdf images
  } else if (tipnum >= 8 & tipnum <= 10){
    omi1 <- 2.5
  } else {
    omi1 <- 2
  }
  return(list(height = hei, omi1 = omi1))

}

#' Wrap a Character String to a Plotting Area.
#'
#' `wrap_string_to_plot()` calculates the optimal `cex` and `width` for a character
#' string to be added as a plot title. It works once a plot device has been called.
#'
#' @details Idea to use [strwrap()] from https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
#' @param string A character vector with the text to be plotted
#' @param max_cex A real number, giving the maximum *cex* (**c**haracter **ex**pansion) for the string to be plotted
#' @param min_cex minimum character expansion to be used on the title
#' @param string_font font type to be used on the title
#' @param max_height A real number, giving the maximum height to be covered by the text
#' @param init_strwrap_width A real number indicating the minimum number of characters to be plotted by line
#' @param min_width A real number, giving the minimum width to be occupied by the string
#' @param max_width A real number, giving the maximum width to be occupied by the string
#' @param whole Boolean, indicating if the whole string should be plotted even if it surpasses the limits established in previous arguments
#' @export
wrap_string_to_plot <- function(string, max_cex = 1, min_cex = 0.5, string_font = 2,
                           max_height = graphics::par("din")[2]-graphics::par("pin")[2]- graphics::par("omi")[1]-graphics::par("mai")[1] - 0.2,
                           init_strwrap_width = 50,  min_width = 0.9*graphics::par("din")[1],
                           max_width = 0.9*graphics::par("din")[1], whole = TRUE){
  # collapse string to a vetor of one element in case it has more elements
  if(max_height <= 0){
    message("string cannot be adjusted because there is not enough space on upper margin of plotting device")
    return(NA)
  }
  wraps <- strwrap(string, width = init_strwrap_width)
  sw <- graphics::strwidth(s = wraps, units = "inches", cex = max_cex, font = string_font)
  sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = max_cex, font = string_font)
  wi <- init_strwrap_width - 1
  string_cex <- max_cex
  # next if when title is too short and fits in one line with regular cex:
  if(sh < max_height & max(sw) < max_width){ #if(length(wraps) ==1)
      return(list(wrapped = paste(wraps, collapse = "\n"), wraps = wraps,
      string_cex = string_cex, strwrap_width = wi, string_font = string_font))
  }
  # next while to find the appropriate cex to fit a long title with a max and min width:
  while(sh > max_height | max(sw) > max_width | max(sw) < min_width) {  #max(sw) < min_width | max(sw) > max_width
    while(max(sw) < min_width) {
      wi <- wi + 1
      wraps <- strwrap(string, width = wi)
      sw <- graphics::strwidth(s = wraps, units = "inches", cex = string_cex, font = string_font)
    }
    sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = string_cex, font = string_font)
    if(sh > max_height){ #length(wraps) > n_lines |
      if(string_cex <= min_cex){
        break
      }
      string_cex <- string_cex - 0.01
      sw <- graphics::strwidth(s = wraps, units = "inches", cex = string_cex, font = string_font)
      sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = string_cex, font = string_font)
    } else {
      break
    }
  }
  if(!whole){
    while(sh > max_height){
      wraps <- wraps[-length(wraps)]
      sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = string_cex, font = string_font)
    }
  }
  return(list(wrapped = paste(wraps, collapse = "\n"), wraps = wraps,
  string_cex = string_cex, strwrap_width = wi, string_font = string_font))
}

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
                           write = "no",
                           folder_name = "phylo_all",
                           file_name = "chronogram",
                           plot_type = "phyloch",
                           max_depth,
                           cex_tiplabels = graphics::par("cex"),
                           cex_axislabel = graphics::par("cex"),
                           cex_axis = graphics::par("cex"),
                           cex_title = graphics::par("cex"),
                           mai1, mai2, mai3, mai4,
                           omi1, omi2, omi3, omi4,
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
  if(missing(mai4))
  mai4 <- unique(unlist(sapply(chronograms, "[", "tip.label")))
  ind <- which.max(nchar(mai4))
  mai4 <- graphics::strwidth(s = mai4[ind],
                             units = "inches",
                             cex = cex_tiplabels,
                             font = 3)
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
  file_prefix <- paste0(folder_name,
                      "/",
                      gsub("\\.png$|\\.pdf$", "", file_name),
                      "_")
  for (i in seq(chronograms)) {
    file_name <- paste0(file_prefix,
                        i,
                        ".",
                        write)
    plot_phylo(chronogram = chronograms[[i]],
               title = names(chronograms)[i],
               time_depth = max_depth,
               plot_type = plot_type,
               cex_tiplabels = cex_tiplabels,
               cex_axislabel = cex_axislabel,
               cex_axis = cex_axis,
               cex_title = cex_title,
               mai4 = mai4,
               write = write,
               file_name = file_name,
               geologic_timescale = NULL,
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

#' @import utils
utils::globalVariables(c("strat2012"))


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
#' @param mai1,mai2,mai3,mai4 A numeric value indicating internal plot margin sizes
#' in inches. `mai4` is particularly important as it indicates the space needed for
#' plotting whole tip labels on the right margin of the plot.
#' @param omi1,omi2,omi3,omi4 A numeric value indicating outter plot margin sizes in inches.
#' @param plot_height,plot_width A numeric value indicating height and width for the plot.
#' @param write A character vector of length 1 indicating the file extension to
#' write the plots to. Options are "pdf" or "png". Anything else will not write a file.
#' @param file_name A character string giving the name and path to write the files to.
#' @param geologic_timescale A dataframe of geochronological limits.
#' @param geologic_unit A character vector used to select geological time units that shall be displayed. When using \code{gradstein04}, \code{"eon"}, \code{"era"}, \code{"period"}, \code{"epoch"}, and \code{"stage"} are available.
#' @param cex_axis A numeric indicating character expansion for the axis. Default to value given by `graphics::par("cex")`.
#' @param cex_title A numeric indicating character expansion for the title. Default to value given by `graphics::par("cex")`.
#' @param cex_axislabel A numeric indicating character expansion for the time axis label. Default to value given by `graphics::par("cex")`.
#' @param pos_title Indicates the line position of the title. Default to 1.
#' @param pos_axis Indicates the line position of the axis. Default to 1.
#' @param center_axislabel A numeric indicating center position of time axis label. Default to 0.5.
#' @param axis_label A character string used to provide information on time units. Defaults
#' to "Time (MYA)" (time in million years ago). If NULL, time units are not added.
#' @param cex_tiplabels A numeric value indicating **c**haracter **ex**pansion (i.e.,
#'  size scaling factor) of tip labels in `chronogram`. Default to value given by `graphics::par("cex")`.
#' @inheritParams phyloch::axisGeo
#' @inheritDotParams ape::plot.phylo
#' @importFrom ape .PlotPhyloEnv
#' @export
# enhance: examples of axis_types!
plot_phylo <- function(chronogram,
                       title = "Chronogram",
                       time_depth = NULL,
                       plot_type = "phyloch",
                       mai1, mai2, mai3, mai4,
                       omi1, omi2, omi3, omi4,
                       plot_height, plot_width,
                       write = "no",
                       file_name = NULL,
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

  if ("png" %in% write) {
    grDevices::png(file = file_name, height = plot_height)
  }
  if ("pdf" %in% write) {
    grDevices::pdf(file = file_name, height = plot_height/72)
  }
  graphics::par(xpd = NA, mai = c(mai1, mai2, mai3, mai4), omi = c(omi1, omi2, omi3, omi4))
  # plot_chronogram.phylo(chronograms[[i]], cex = 1.5, edge.width = 2, label.offset = 0.5,
    # x.lim = c(0, max_depth), root.edge = TRUE, root.edge.color = "white")
  # graphics::par(xpd = FALSE)
  if (plot_type %in% c("ape", "phyloch")) {
    ape::plot.phylo(chronogram,
                    cex = cex_tiplabels, #edge.width = 2,
                    label.offset = 0.5,
                    x.lim = c(0, time_depth),
                    root.edge = TRUE,
                    plot = TRUE, ...)  #
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
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
  if (!is.null(units)) {
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
  if (any(c("png", "pdf") %in% write)) {
    grDevices::dev.off()
  }
}
# chronogram <- plant_bold_otol_tree
# plot_phylo_gg <- function(chronogram, title = "Tree", time_depth = NULL, plot_type = 1,
# cex = graphics::par("cex"), mai4 = NULL, write = "nothing", file_name = NULL, GTS = utils::getAnywhere("strat2012")){
#   max_age <- max(ape::branching.times(chronogram))
#   age_lim <- max_age*1.2
#   grDevices::pdf("test.pdf")
#   p <- ggtree::ggtree(chronogram) + ggtree::geom_tiplab()  + #ggplot2::xlim(age_lim*0.1,-age_lim) +
#   ggplot2::coord_cartesian(xlim = c(age_lim*0.5,-age_lim), ylim = c(-1, ape::Ntip(chronogram)), expand = FALSE) +
#   ggplot2::scale_x_continuous(breaks=seq(-age_lim,0,100), labels=abs(seq(-age_lim,0,100))) +
#   ggtree::theme_tree2()
#   p <- ggtree::revts(p)
#   deeptime::gggeo_scale(p, neg = TRUE)
#   print(p)
#   grDevices::dev.off()
# }
# .PlotPhyloEnv <- new.env()
