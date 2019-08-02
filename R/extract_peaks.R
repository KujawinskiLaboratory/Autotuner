#' @title extract_peaks
#'
#' @description This function is designed to extract peaks observed within the
#' TIC from each sample, and return their indicies for further processing.
#'
#' @param Autotuner An Autotuner objected containing sample specific raw
#' data.
#' @param returned_peaks A scalar number of peaks to return for visual
#' inspection. Five is the minimum possible value.
#' is the standard.
#' @param signals A list containing traces and locations where signals are
#' detected across all samples being checked by the algorithm.
#'
#' @return peak_table_list - a list of data.frame tables containing information
#' on where each where peaks are located within each sample.
extract_peaks <- function(Autotuner,
                          returned_peaks = 10,
                          signals) {


    factorCol <- Autotuner@factorCol
    metadata <- Autotuner@metadata

    assertthat::assert_that(is.character(factorCol),
                          msg = "Make sure factorCol is a string object.")

    if(!(factorCol %in% colnames(metadata))) {
        stop("factorCol is not within the Autotuner metadata.")
    }


    sample_names <- paste(unlist(metadata[,factorCol]),
                          1:nrow(metadata))

    peak_table_list <- list()

    # Running Algorithm to extract peaks from each sample ---------------------
    for(index in seq_along(sample_names)) { #making peak table for each sample

        ## identifying regions within TIC where peaks were identified

        #### 2019-07-02: extending the window to make sure all peaks have
        #### atleast two scans
        # <- which(signals[[index]]$signals == 1)
        peaks <- rle(signals[[index]]$signals)
        counter <- 1
        peakGroups <- list()
        for(rleIndex in seq_along(peaks$lengths)) {

            curValue <- peaks$values[rleIndex]
            if(curValue == 1) {

                peakGroups[[counter]] <- data.frame(index = (startPoint + 1),
                                                    length =
                                                        peaks$lengths[rleIndex])
                startPoint <- startPoint + peaks$lengths[rleIndex]
                counter <- counter + 1

            } else {

                if(rleIndex == 1) {
                    startPoint <- peaks$lengths[rleIndex]
                } else {
                    startPoint <- startPoint + peaks$lengths[rleIndex]
                }


            }

        }
        peakGroups <- Reduce(peakGroups, f = rbind)
        signals[[index]]$signals[peakGroups$index[peakGroups$length
                                                  == 1] + 1] <- 1

        #peaks <- Autotuner@time[[index]][signals[[index]]$signals == 1]

        #peaks <- peaks[!is.na(peaks)]

        ## getting and extracting peaks
        #peaks
        findPeaks <- rle(signals[[index]]$signals == 1)
        counter <- 1
        peakGroups <- list()
        for(rleIndex in seq_along(findPeaks$lengths)) {

            curValue <- findPeaks$values[rleIndex]
            if(curValue) {

                start <- (startPoint + 1)
                startPoint <- startPoint + findPeaks$lengths[rleIndex]
                end <- startPoint
                peakGroups[[counter]] <- data.frame(start,
                                                    end,
                                                    length = end - start + 1)
                counter <- counter + 1

            } else {

                if(rleIndex == 1) {
                    startPoint <- findPeaks$lengths[rleIndex]
                } else {
                    startPoint <- startPoint + findPeaks$lengths[rleIndex]
                }


            }

        }

        peakGroups <- Reduce(rbind, peakGroups)

        ## generating pseudo retention time correction threshold
        ## based on the distribution of TIC peaks within the sample

        ### 2019-07-02 MOVING AWAY FROM THRESHOLDS SINCE THEY DON'T SEEM TO
        ### APPLY ANY LONGER. SEEMS LIKE I USED THEM IN THE PAST TO ASIGN THINGS
        ### THAT MAY NOT HAVE HAD ADJACENT SCANS. REMOVING IT WOUDL BE BETTER.
        #threshold <- estimate_threshold(distance_vector = peak_dist)


        peakGroups
        if(nrow(peakGroups) < returned_peaks) {
            returned_peaks <- nrow(peakGroups)
        }
        peakGroups <- peakGroups[order(peakGroups$length, decreasing = TRUE),]

        peak_times <- list()
        for(j in 1:nrow(peakGroups)) {

            peak_times[[j]] <- Autotuner@time[[index]][
                peakGroups$start[j]:peakGroups$end[j]]

        }


        max_peak_length <- max(vapply(X = peak_times, FUN = length,
                                      FUN.VALUE = numeric(1)))
        peak_table <- data.frame(matrix(nrow = max_peak_length,
                                            ncol = returned_peaks+1))
        peak_table[,1] <- 1:max_peak_length
        colnames(peak_table) <- c("peakLenth", paste("peak", 1:returned_peaks))

        for(column in 2:ncol(peak_table)) {
            peak <- peak_times[[column -1]]
            peak_table[c(1:length(peak)),column] <- peak
        }

        peak_table$peakLenth <- NULL
        peak_table_list[[index]] <- peak_table


    }

    names(peak_table_list) <- sample_names
    rm(index)
    return(peak_table_list)
}

