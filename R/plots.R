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

#' subset trees for plotting densitree plots and phylo_all plots
#' @param trees A list of trees as multiPhylo or as a plain list object.
#' @param include Boolean or numeric vector. Default to TRUE, keep all chronograms
#' in trees. If FALSE, exclude chronograms with only two tips. If numeric, it is used
#' as indices to subset trees object.
subset_trees <- function(trees, include = TRUE){
  if(is.numeric(include)){  #if it is numeric
    include <- round(include)
    include <- include[which(include <= length(trees))]
    include <- sort(unique(include))
    trees <- trees[include]
  }
  if(is.logical(include) & !is.na(include)){
      if (!include){
      trees <- trees[sapply(trees, function (x) ape::Ntip(x)) > 2]
    }
  }
  trees
}

#' get a densiTree plot from a set of opentree_chronograms
#' if densiTree plot function throws an error, it chooses the tree with the most tips as consensus (using get_biggest_phylo)
#' we found that densiTree errors commonly from failing to do a consensus tree.
#' @param trees A list of trees as multiPhylo or as a plain list object.
#' @inheritParams subset_trees
#' @inheritDotParams phangorn::densiTree -x -consensus
#' @export
plot_densitree <- function(trees, include = TRUE, ...){
  trees <- subset_trees(trees, include = include)
  # for densitree plot it does not really matter if trees have the same length
  # max_depth <- round(max(sapply(trees, function(x) max(ape::branching.times(x)))) + 5, digits = -1)
  # # max_depth <- round(max(summarize_datelife_result(datelife_result = get.filtered.results(),
  # #             summary_format = "mrca")) + 5)  # in case we used it for a datelife object
  #
  # # plot all trees from the same depth:
  # trees <- lapply(trees, function(x) {
  #  x$root.edge <- max_depth - max(ape::branching.times(x))
  #  x
  # })
  class(trees) <- "multiPhylo"
  # if we use biggest phylo as consensus for all, some data set are not plotted correctly
  # biggest_phylo <- datelife::get_biggest_multiphylo(trees = trees)
  # try(phangorn::densiTree(x = trees, consensus = biggest_phylo, ...))
  tryCatch(phangorn::densiTree(x = trees, ...),
  error = function(e) {
    biggest_phylo <- datelife::get_biggest_multiphylo(trees = trees)
    try(phangorn::densiTree(x = trees, consensus = biggest_phylo, ...))
  })
}

#' get the outer margin of a graphics device from a number of tips
#' @param tree phylo object to be plotted
phylo_height_omi <- function(tree){
  tipnum <- ape::Ntip(tree)
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

#' wrap a character string to a plotting area. It gives the optimal cex and width. It works once plot device has been called
#' idea to use strwrap from https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
#'
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
#' @param trees A list of chronograms as `multiPhylo` or as a plain `list` object.
#' @inheritParams subset_trees
#' @inheritParams plot_phylo
#' @param individually Boolean indicating if trees should be plotted one by one
#' or appended to the same file.
#' @inheritDotParams ape::plot.phylo
#' @param folder_name A character string indicating the name of the folder to write
#' the chronogram plot files to. Default to "phylo_all". Only relevant if `write = "png"` or `write = "pdf"`.
#' @param file_name A character string indicating a file name to write the chronogram
#' plots to. Relevant if `write = "png"` or `write = "pdf"`.
#' @param max_depth A numeric vector of length 1, indicating the upper limit of
#'   the time scale on the x axis to be used on all plots. If none is provided, it is estimated
#'   by getting the tree depth of the oldest chronogram, adding 5 and rounding to the
#'   closest 10. See details for more.
#' @details
#' Currently, max_depth is obtained by default with `round(max(sapply(trees, function(x) max(ape::branching.times(x)))) + 5, digits = -1)`.
#' @export
plot_phylo_all <- function(trees,
                           cex = graphics::par("cex"),
                           include = TRUE,
                           individually = TRUE,
                           write = "no",
                           folder_name = "phylo_all",
                           file_name = "chronogram",
                           plot_type = "phyloch",
                           max_depth, ...) {
  trees <- subset_trees(trees, include = include)
  # in case there is just one tree in trees
  if (any("tip.label" %in% names(trees))) {
    trees <- list(trees)
  }
  if (missing(max_depth)) {
    max_depth <- round(max(sapply(trees,
                                  function(x) max(ape::branching.times(x)))) + 5, digits = -1)
    # if(isTRUE(all.equal(round(sapply(trees,
    # function(x) max(ape::branching.times(x))), digits = 3))))
    # max_depth <- round(max(summarize_datelife_result(datelife_result = get.filtered.results(),
    #             summary_format = "mrca")) + 5)  # in case we used it for a datelife object

  }
  if (!is.numeric(max_depth)) {
    stop("'max_depth' argument must be 'numeric'. It currently is ", max_depth)
  }
  # plot all trees from the same depth:
  trees <- lapply(trees, function(x) {
   x$root.edge <- max_depth - max(ape::branching.times(x))
   x
  })
  mai4 <- unique(unlist(sapply(trees, "[", "tip.label")))
  ind <- which.max(nchar(mai4))
  mai4 <- graphics::strwidth(s = mai4[ind],
                             units = "inches",
                             cex = cex,
                             font = 3)
  # if(any(lapply(trees, ape::Ntip) > 3))
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
  for (i in seq(trees)) {
    file_name <- paste0(file_prefix,
                        i,
                        ".",
                        write)
    plot_phylo(trees[[i]],
               names(trees)[i],
               time_depth = max_depth,
               plot_type = plot_type,
               cex = cex,
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
#' of a given tree with a time axis representing geologic time.
#'
#' @param tree A chronogram either as a newick character string or as a `phylo` object.
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
#' @param mai4 A numeric vector of length one indicating the space needed for
#'   plotting whole tip labels (right margin of the plot).
#' @param write A character vector of length 1 indicating the file extension to
#' write the plots to. Options are "pdf" or "png". Anything else will not write a file.
#' @inheritParams phyloch::axisGeo
#' @param file_name A character string giving the name and path to write the files to.
#' @param geologic_timescale A dataframe of geochronological limits.
#' @inheritDotParams ape::plot.phylo
#' @export
# enhance: examples of axis_types!
plot_phylo <- function(tree,
                       title = "Tree",
                       time_depth = NULL,
                       plot_type = "phyloch",
                       cex = graphics::par("cex"),
                       mai4 = NULL,
                       write = "no",
                       file_name = NULL,
                       geologic_timescale = "strat2012",
                       unit = "period",
                       ...){
  #
  if (!inherits(tree, "phylo")) {
      if (inherits(tree, "multiPhylo")) {
        message("'tree' is a 'multiPhylo' object. Only first tree will be plotted.")
        tree <- tree[[1]]
      } else {
        message("'tree' must be a 'phylo' or 'multiPhylo' object.")
        return(NA)
      }
  }
  if (is.null(tree$edge.length)) {
    message("Tree has no edge lengths; a tree with a geologic time axis can not be plotted.")
    return(NA)
  }
  phylo_length <- max(ape::branching.times(tree))
  if (is.null(time_depth)) {
    if (is.null(tree$root.edge)) {
      time_depth <- round(phylo_length*1.2, digits = -1)
    } else {
      time_depth <- round(phylo_length + tree$root.edge, digits = -1)
    }
  }
  match.arg(arg = plot_type, choices = c("phyloch", "strap", "phytools", "ape"))
  if (is.null(geologic_timescale) | "strat2012" %in% geologic_timescale) {
    utils::data("strat2012", package = "phyloch")
    geologic_timescale <- strat2012
  }
  if (is.null(mai4)) {
    ind <- which.max(nchar(tree$tip.label))
    mai4 <- graphics::strwidth(s = tree$tip.label[ind],
                               units = "inches",
                               cex = cex,
                               font = 3)
  }
  pho <- phylo_height_omi(tree = tree)
  if ("png" %in% write) {
    grDevices::png(file = file_name, height = pho$height)
  }
  if ("pdf" %in% write) {
    grDevices::pdf(file = file_name, height = pho$height/72)
  }
  graphics::par(xpd = NA, mai = c(0, 0, 0, mai4), omi = c(pho$omi1, 0, 1, 0))
  # plot_chronogram.phylo(trees[[i]], cex = 1.5, edge.width = 2, label.offset = 0.5,
    # x.lim = c(0, max_depth), root.edge = TRUE, root.edge.color = "white")
  # graphics::par(xpd = FALSE)
  if (plot_type %in% c("ape", "phyloch")) {
    ape::plot.phylo(tree,
                    cex = cex, #edge.width = 2,
                    label.offset = 0.5,
                    x.lim = c(0, time_depth),
                    root.edge = TRUE,
                    plot = TRUE, ...)  #
    if ("ape" %in% plot_type) {
      ape::axisPhylo()
    } else { # if ("phyloch" %in% plot_type) {
      axisGeo(GTS = geologic_timescale,
              unit = unit,
              col = c("gray80", "white"),
              gridcol = c("gray80", "white"),
              cex = 0.5,
              gridty = "twodash")
    }
    graphics::mtext("Time (MYA)",
                    cex = cex,
                    side = 1,
                    font = 2,
                    line = (pho$omi1-0.2)/0.2,
                    outer = TRUE,
                    at = 0.5)
  }
  if ("phytools" %in% plot_type) {
    # TODO
    message("Plotting a geologic time axis with phytools is not supported yet.")
  }
  if ("strap" %in% plot_type) {
    tree$root.time <- phylo_length
    strap::geoscalePhylo(tree = tree,
                         x.lim = c(0, phylo_length),
                         cex.tip = 0.7,
                         cex.ts = 0.7,
                         cex.age = 0.7,
                         width = 4,
                         tick.scale = 15,
    # creates boxes with the last unit in argument "unit":
                         boxes = unit[length(unit)],
                         erotate = 90,
                         quat.rm = TRUE,
                         units = unit)
    graphics::mtext("Time (MYA)",
                    cex = cex,
                    side = 1,
                    font = 2,
                    line = (pho$omi1-0.2)/0.2,
                    outer = FALSE,
                    at = 1)
  }
  # add a title to the plot
  if (!is.null(title)) {
    titlei <- wrap_string_to_plot(string = title, max_cex = 1, whole = FALSE)
    graphics::mtext(text = titlei$wrapped, outer = TRUE,
      cex = titlei$string_cex, font = titlei$string_font, line = 1)
  }
  if (any(c("png", "pdf") %in% write)) {
    grDevices::dev.off()
  }
}
# tree <- plant_bold_otol_tree
# plot_phylo_gg <- function(tree, title = "Tree", time_depth = NULL, plot_type = 1,
# cex = graphics::par("cex"), mai4 = NULL, write = "nothing", file_name = NULL, GTS = utils::getAnywhere("strat2012")){
#   max_age <- max(ape::branching.times(tree))
#   age_lim <- max_age*1.2
#   grDevices::pdf("test.pdf")
#   p <- ggtree::ggtree(tree) + ggtree::geom_tiplab()  + #ggplot2::xlim(age_lim*0.1,-age_lim) +
#   ggplot2::coord_cartesian(xlim = c(age_lim*0.5,-age_lim), ylim = c(-1, ape::Ntip(tree)), expand = FALSE) +
#   ggplot2::scale_x_continuous(breaks=seq(-age_lim,0,100), labels=abs(seq(-age_lim,0,100))) +
#   ggtree::theme_tree2()
#   p <- ggtree::revts(p)
#   deeptime::gggeo_scale(p, neg = TRUE)
#   print(p)
#   grDevices::dev.off()
# }
# .PlotPhyloEnv <- new.env()
