#' @title plot_signals
#'
#' @description this funciton plots the peak identified within chromatography.
#'
#' @param Autotuner Autotuner object created following Create_Autotuner()
#' initialization function.
#' @param threshold User input scalar value for the number of standard
#' deviations required to consider a peak to be significant.
#' @param sample_index which of all of the samples should the user plot.
#' Entered in as a numerical index value with length 1.
#' @param signals A vector containing information on where signals are
#' located between samples.
#'
#'
#' @details This function plots the chromatography and the matching sliding
#' window signal processing results for each sample. Signal processing
#' functions will be the same color as the chromatography spectra, just a
#' lighter shade and a different type of line. The chromatography will be a
#' solid line, the signal will be a dashed line (lty = 2) with a .5 alpha,
#' and the thersholds will be dotted lines (lty = 3) with alpha values of 3.
#'
#' @return Plots signal
#'
#' @examples
#' Autotuner <- readRDS(system.file("extdata/Autotuner.rds",
#' package="Autotuner"))
#' lag <- 25
#' threshold <- 3.1
#' influence <- 0.1
#'
#' signals <- lapply(getAutoIntensity(Autotuner),
#' ThresholdingAlgo, lag, threshold, influence)
#'
#' plot_signals(Autotuner, threshold, sample_index = seq_len(3), signals = signals)
#'
#' @export
plot_signals <- function(Autotuner, threshold, sample_index, signals) {

    factorCol <- getAutoFactorCol(Autotuner)
    metadata <- getAutoMetadata(Autotuner)


    # checking for correct input ----------------------------------------------
    if(length(sample_index) > 8) {
        sample_index <- sample_index[seq_len(8)]
    }

    if(!(factorCol %in% colnames(metadata))) {
        stop("factorCol is not within the MSest metadata.")
    }

    if(length(signals) == 0) {
        stop("Signals for the samples have not been computed.")
    }



    # generating plots --------------------------------------------------------
    sample_type <- metadata[sample_index,
                            grep(factorCol,colnames(metadata))]
    par(mfcol = c(2,1)
        #oma = c(0,2,0,0) + 0.1,
        #mar = c(1,2,1.5,1)
        )

    if(length(sample_index) < 3) {
        cols <- c("red", "blue")
    } else {
        cols <- RColorBrewer::brewer.pal(length(sample_index), name = "Dark2")
    }


    # ploting chromatogram ----------------------------------------------------
    par(mai = c(.2,.8,.3,.1))
    for(i in seq_along(sample_index)) { # plotting chromatograms together

        # removed secondary plots - should do one plot per sample groups
        index <- sample_index[i]
        curSignal <- signals[[i]]
        curCol <- cols[i]

        ave_filter <- curSignal$avgFilter
        signal_subset <- seq_along(ave_filter)


        # plotting TICs ------------------------------------------------------
        if(i == 1) {
            plot(Autotuner@time[[index]],
                 Autotuner@intensity[[index]],
                 type="l",
                 ylab="Intensity",
                 xlab="",
                 main = paste("TIC Peak Detection Across Samples"),
                 xaxt='n',
                 col = scales::alpha(curCol, 1))
        } else {
            lines(Autotuner@time[[index]],
                  Autotuner@intensity[[index]],
                  type="l",
                  col = scales::alpha(curCol, 1))

        }


        # ploting signal processing function --------------------------------
        lines(Autotuner@time[[index]][signal_subset],
              ave_filter,
              type="l",
              col=scales::alpha(curCol, 0.5),
              lty= 1,
              lwd = 3)
        lines(Autotuner@time[[index]][signal_subset],
              ave_filter +
                threshold *
                curSignal$stdFilter,
              type="l",
              col=scales::alpha(curCol, 0.5),
              lty=3,
              lwd = 3)

    }


    # plotting signal plots below chromatogram -----------------------------
    allSignals <- list()
    for(i in seq_along(sample_index)) {

        # idea - align signals over time and highlight regions over lap in a
        # single plot
        plotSignals <- signals[[i]]$signals
        signalDf <- data.frame(time = Autotuner@time[[i]][
            seq_along(plotSignals)],
                               signals = plotSignals)
        rm(plotSignals)

        signalDf$signals[signalDf$signals == -1] <- 0
        colnames(signalDf) <- paste(colnames(signalDf), i, sep = "_")
        allSignals[[i]] <- signalDf

    }

    minRow <-  min(vapply(X = allSignals,
                          FUN = nrow,
                          FUN.VALUE = numeric(1)))
    allSignals <- lapply(allSignals, function(x) {x[seq_len(minRow),]})

    if(length(sample_index) == 1) {
        combinedSignals <- allSignals[[1]]
        overlappingSignals <- combinedSignals$signals_1
    } else {
        combinedSignals <- Reduce(cbind, allSignals)
        overlappingSignals <- combinedSignals[,grep("signal",
                                                    colnames(combinedSignals))]
        overlappingSignals <- rowSums(overlappingSignals) > 1

    }
    par(mai = c(.7,.8,.2,.1))
    plot(combinedSignals$time_1,
         overlappingSignals,
         type="S",
         ylab = "Peak Overlap",
         xlab="",
         ylim=c(0,1.2),
         lwd=1,
         las=1)
    graphics::title(xlab="Time (s)", line=2.1, cex.lab=1.2)
    legend("topright", fill = cols,
           legend = sample_type,
           title = "Samples Visualized",
           horiz=TRUE, cex=0.7)


}
