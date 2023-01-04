#' Plot Spatial Point Process
#'
#' This function is used to plot a spatial point process.
#'
#' @param data Data.frame with x,y, and agent type (in that order)
#' @param colorGuide (Optional) String for guides(color=) in ggplot2. Usually
#'     NULL or 'none' is sufficient. Default is NULL.
#' @param ptSize (Optional) Numeric indicating point size. Default is 1.
#' @param xlim (Optional) Two value numeric vector indicating the size of the
#'     region in the x-direction. Default is c(min(x),max(x)).
#' @param ylim (Optional) Two value numeric vector indicating the size of the
#'     region in the y-direction. Default is c(min(y),max(y)).
#' @param dropAxes (Optional) Boolean indicating if the x, y axis title and
#'     labels should be dropped. Default is FALSE.
#' @param colors (Optional) Vector of colors for the points. Default is NULL, or
#'     ggplot2 selected colors.
#'
#' @return ggplot2 plot of the spatial point process.
#' @export
#'
#' @examples
#' dat <- simulatePP(cellVarData=
#'                        data.frame('stage'=c(0,1),
#'                                   'A'=c(0,0),
#'                                   'B'=c(1/50,1/50)),
#'                    cellKappaData=data.frame(
#'                                   'cell'=c('A','B'),
#'                                   'clusterCell'=c(NA,'A'),
#'                                   'kappa'=c(20,5)),
#'                    peoplePerStage=100,
#'                    imagesPerPerson=1,
#'                    reduceEdge=0.025,
#'                    silent=F )
#' plotPP(dat[dat$Image==1,c('x','y','cellType')])
plotPP <- function(data, colorGuide = NULL,ptSize=1,
                   xlim=c(min(data[,1]),max(data[,1])),
                   ylim=c(min(data[,2]),max(data[,2])),
                   dropAxes=F,
                   colors=NULL){
  ## This is used to plot point process data

  retPlot <- ggplot2::ggplot() +
    ggplot2::geom_point(mapping=ggplot2::aes(x=data[,1],y=data[,2],
                                             col=data[,3]),
               size=ptSize) +
    ggplot2::theme_bw() +
    ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim) +
    ggplot2::guides(color=colorGuide)
  if(dropAxes)
    retPlot <- retPlot +
      ggplot2::theme(axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank())
  if(!is.null(colors))
    retPlot <- retPlot +
      ggplot2::scale_color_manual(values = colors)

  retPlot
}
