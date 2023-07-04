#' Plot Spatial Point Process
#'
#' This function is used to plot a spatial point process.
#'
#' @param data Data.frame with x,y, and agent type (in that order)
#' @param colorGuide (Optional) String for guides(color=) in ggplot2. Usually
#'     NULL or 'none' is sufficient, but ggplot2::guide_legend() can also be
#'     used for more custom results. Default is NULL.
#' @param ptSize (Optional) Numeric indicating point size. Default is 1.
#' @param xlim (Optional) Two value numeric vector indicating the size of the
#'     region in the x-direction. Default is c(min(x),max(x)).
#' @param ylim (Optional) Two value numeric vector indicating the size of the
#'     region in the y-direction. Default is c(min(y),max(y)).
#' @param dropAxes (Optional) Boolean indicating if the x, y axis title and
#'     labels should be dropped. Default is FALSE.
#' @param layerBasedOnFrequency (Optional) Boolean indicating if the data should be
#'     layer based on the number of cells of the type. Default is TRUE.
#' @param colors (Optional) Vector of colors for the points. Default is NULL, or
#'     ggplot2 selected colors.
#'
#' @return ggplot2 plot of the spatial point process.
#' @export
#'
#' @examples
#' plotPP(TNBC_pheno[TNBC_pheno$Person == 1, c("cellx", "celly", "Phenotype")])
#' plotPP(diabetes[diabetes$Image == "E37", c("x", "y", "cellType")])
plotPP <- function(data, colorGuide = NULL, ptSize = 1,
                   xlim = c(min(data[, 1]), max(data[, 1])),
                   ylim = c(min(data[, 2]), max(data[, 2])),
                   dropAxes = FALSE, layerBasedOnFrequency = TRUE,
                   colors = NULL) {
  # Sort so most populous cells are at the bottom
  if (layerBasedOnFrequency && length(unique(data[, 3])) > 1) {
    cells_order <- as.data.frame(table(data[, 3])[order(-table(data[, 3]))])$Var1

    idxs <- sapply(cells_order, function(x, data1) {
      which(data1[, 3] == x)
    }, data1 = data)
    data <- data[unlist(idxs), ]
    rownames(data) <- NULL
  }

  # Plot
  retPlot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = data[, 1], y = data[, 2],
        col = data[, 3]
      ),
      size = ptSize
    ) +
    ggplot2::theme_bw() +
    ggplot2::xlim(xlim) +
    ggplot2::xlab(NULL) +
    ggplot2::ylim(ylim) +
    ggplot2::ylab(NULL) +
    ggplot2::guides(color = colorGuide)
  if (dropAxes) {
    retPlot <- retPlot +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank()
      )
  }
  if (!is.null(colors)) {
    retPlot <- retPlot +
      ggplot2::scale_color_manual(values = colors)
  }

  retPlot
}


#' Compare K Functions Between Outcomes
#'
#' @param data Data.frame with the columns r, K, Unit, and Outcome. r indicates
#'   the radius of checked K function, K indicates the K function value, Unit
#'   specifies the unique K function, and Outcome indicates the unit outcome.
#' @param inc.legend (Optional) Boolean indicating if the legend should be given.
#'   This will also include numbers to indicate if any K functions are missing.
#'   The default is TRUE.
#' @param inc.noise (Optional) Boolean indicating if a gray, dashed line should
#'   be included to show what total noise would be like. The default is FALSE.
#'
#' @return ggplot2 object showing the K function with a superimposed average
#' @export
#'
#' @examples
#' # Example 1
#' tmp <- getKFunction(TNBC_pheno[TNBC_pheno$Class == 0, -1],
#'   agents = c("Tumor", "Tumor"),
#'   unit = "Person",
#'   rCheckVals = seq(0, 50, 1)
#' )
#' tmp1 <- getKFunction(TNBC_pheno[TNBC_pheno$Class == 1, -1],
#'   agents = c("Tumor", "Tumor"),
#'   unit = "Person",
#'   rCheckVals = seq(0, 50, 1)
#' )
#' tmp_1 <- tidyr::pivot_longer(data = tmp, cols = K1:K18)
#' tmp1_1 <- tidyr::pivot_longer(data = tmp1, cols = K1:K15)
#'
#' data_plot <- rbind(
#'   data.frame(
#'     "r" = tmp_1$r,
#'     "K" = tmp_1$value,
#'     "Unit" = tmp_1$name,
#'     "Outcome" = "0"
#'   ),
#'   data.frame(
#'     "r" = tmp1_1$r,
#'     "K" = tmp1_1$value,
#'     "Unit" = paste0(tmp1_1$name, "_1"),
#'     "Outcome" = "1"
#'   )
#' )
#'
#' plot_K_functions(data_plot)
#'
#' # Example 2
#' tmp <- getKFunction(TNBC_pheno[TNBC_pheno$Class == 0, -1],
#'   agents = c("Tumor", "B"), unit = "Person",
#'   rCheckVals = seq(0, 50, 1)
#' )
#' tmp1 <- getKFunction(TNBC_pheno[TNBC_pheno$Class == 1, -1],
#'   agents = c("Tumor", "B"), unit = "Person",
#'   rCheckVals = seq(0, 50, 1)
#' )
#'
#' tmp_1 <- tidyr::pivot_longer(data = tmp, cols = K1:K18)
#' tmp1_1 <- tidyr::pivot_longer(data = tmp1, cols = K1:K15)
#'
#' data_plot <- rbind(
#'   data.frame(
#'     "r" = tmp_1$r,
#'     "K" = tmp_1$value,
#'     "Unit" = tmp_1$name,
#'     "Outcome" = "0"
#'   ),
#'   data.frame(
#'     "r" = tmp1_1$r,
#'     "K" = tmp1_1$value,
#'     "Unit" = paste0(tmp1_1$name, "_1"),
#'     "Outcome" = "1"
#'   )
#' )
#'
#' plot_K_functions(data_plot)
plot_K_functions <- function(data, inc.legend = TRUE, inc.noise = FALSE) {
  # Add this to remove notes when building package
  r <- K <- Unit <- Outcome <- Value <- NULL

  if (inc.legend) {
    info <- unique(data[, c("Unit", "Outcome")])
    info$Missing <- NA
    for (i in 1:nrow(info)) {
      info[i, 3] <- nrow(data[data$Unit == info$Unit[i] &
        !stats::complete.cases(data), ]) > 0
    }
    # info$Complete <- data[!stats::complete.cases(data),]
    # info_missing <- unique(data[!stats::complete.cases(data),c('Unit','Outcome')])

    info_labels <- data.frame(
      "Outcome" = unique(data$Outcome),
      "Units" = NA,
      "Missing" = NA,
      "Label" = NA
    )
    for (i in 1:nrow(info_labels)) {
      info_labels[i, 2] <- nrow(info[info$Outcome == info_labels[i, 1], ])
      info_labels[i, 3] <- nrow(info[info$Outcome == info_labels[i, 1] &
        info$Missing, ])
      # Only label Outcomes with numbers if at least one is missing in whole set
      if (sum(info$Missing) != 0) {
        info_labels[i, 4] <- paste0(
          info_labels[i, 1], " (",
          info_labels[i, 2] - info_labels[i, 3], "/",
          info_labels[i, 2], ")"
        )
      } else {
        info_labels[i, 4] <- info_labels[i, 1]
      }
    }
  }

  # Build Averages
  data_wide <- tidyr::pivot_wider(data,
    names_from = "Unit",
    values_from = "K"
  )
  data_avg <- data.frame("r" = unique(data_wide$r))

  for (i in 1:nrow(info_labels)) {
    data_avg[[info_labels[i, "Outcome"]]] <-
      rowMeans(data_wide[data_wide$Outcome == info_labels[i, "Outcome"], -c(1:2)],
        na.rm = TRUE
      )
  }
  data_avg <- tidyr::pivot_longer(data_avg, cols = -r)
  colnames(data_avg) <- c("r", "Outcome", "Value")

  return_plot <-
    ggplot2::ggplot(
      data = stats::na.omit(data),
      ggplot2::aes(
        x = r, y = K,
        group = Unit, color = Outcome
      )
    ) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = r, y = Value,
        group = Outcome, color = Outcome
      ),
      data = data_avg, linewidth = 1.25
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  if (inc.legend) {
    return_plot <- return_plot +
      ggplot2::theme(legend.position = "right") +
      ggplot2::scale_colour_discrete(
        labels = info_labels$Label,
        name = paste0(
          "Outcome\n( ",
          sum(info_labels$Units) - sum(sum(info_labels$Missing)), " / ",
          sum(info_labels$Units), " )"
        )
      )
  }

  if (inc.noise) {
    return_plot <- return_plot +
      ggplot2::geom_line(ggplot2::aes(x = r, y = pi * r^2),
        data = data_avg, linewidth = 1.25,
        color = "gray", linetype = "dashed"
      )
  }

  return_plot
}
