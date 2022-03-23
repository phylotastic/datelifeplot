#' A data frame containing the stratigraphic chart issued in 2012 by the International Commission on Stratigraphy.
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
