#' Axis for time-calibrated phylogenies
#'
#' This function adds a time scale on the side of a phylogeny plot.
#' @param side A numeric value specifying the side where the axis is plotted:
#' 1: below, 2: left, 3: above, 4: right.
#' @param unit A character string used to provide information on time units; defaults to "Ma" (million years ago).
#' @param fact A real number used to scale the axis.
#' @inheritDotParams graphics::axis
#' @author Christoph Heibl
#' @export
axisChrono <- function (side = 1, unit = NULL, fact = 1, ...) {
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (lastPP$type %in% c("phylogram", "cladogram")) {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            x <- pretty(lastPP$xx)
            if (lastPP$direction == "rightwards")
                maxi <- max(lastPP$xx)
            else {
                maxi <- min(lastPP$xx)
                x <- -x
            }
        }
        else {
            x <- pretty(lastPP$yy)
            if (lastPP$direction == "upwards")
                maxi <- max(lastPP$yy)
            else {
                maxi <- min(lastPP$yy)
                x <- -x
            }
        }
    }
    graphics::axis(side = side, at = c(maxi - x), labels =
        abs(x * fact), ...)
    if(!is.null(unit)) {
      graphics::mtext(text = unit, side = side, at = 1.07 * maxi, line = 1, ...)
    }
}
