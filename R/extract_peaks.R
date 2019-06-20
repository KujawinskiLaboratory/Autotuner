#' @title extract_peaks
#'
#' @description This function is designed to extract peaks observed within the
#' TIC from each sample, and return their indicies for further processing.
#'
#' @param Autotuner - An Autotuner objected containing sample specific raw
#' data.
#' @param returned_peaks - A scalar number of peaks to return for visual
#' inspection. Five is the minimum possible value.
#' is the standard.
#' @param signals - List containing traces and locations where signals are
#' detected across all samples being checked by the algorithm.
#'
#' @return peak_table_list - a list of data.frame tables containing information
#' on where each where peaks are located within each sample.
#'
#' @export
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
        peaks <- Autotuner@time[[index]][signals[[index]]$signals == 1]

        peaks <- peaks[!is.na(peaks)]
        ## generating the distance between peaks
        peak_dist <- diff(peaks)
        # rle is done on the difference vector - missing first value

        ## generating pseudo retention time correction threshold
        ## based on the distribution of TIC peaks within the sample
        threshold <- estimate_threshold(distance_vector = peak_dist)


        if(threshold > 0) {

            ## checking for peaks that may come from the same peak.
            run_length <- rle(peak_dist < threshold)
            names(run_length$lengths)[length(run_length$lengths)] <-
                names(run_length$values)[length(run_length$lengths)]

            ## checking to make sure there are enough peaks to return
            passed_threshold <- which(run_length$values == T)
            if(length(passed_threshold) < returned_peaks) {
                returned_peaks <- length(passed_threshold)
            }

            # getting peaks that appear to come from single TIC peak
            top_runs <- sort(run_length$lengths[passed_threshold],decreasing = T)
            top_runs <- match(names(top_runs),
                              names(run_length$lengths))[1:returned_peaks]

            peak_points <- list()

            for(i in 1:length(top_runs)) { # extracting peak info between files

                run <- top_runs[i]

                # correcting for diff
                end <- which(names(peaks) %in% names(run_length$values[run]))

                if(length(run) > 0 && run == 1) {
                    start <- 1
                } else {
                    start <- which(names(peaks) %in%
                                       names(run_length$lengths[run-1]))
                }
                peak_points[[i]] <- peaks[start:end]
            }

            ## store the data
            max_peak_length <- max(sapply(peak_points, length))

            peak_table <- data.frame(matrix(nrow = max_peak_length,
                                            ncol = returned_peaks+1))
            peak_table[,1] <- 1:max_peak_length
            colnames(peak_table) <- c("peakLenth", paste("peak", 1:returned_peaks))

            for(column in 2:ncol(peak_table)) {
                peak <- peak_points[[column -1]]
                peak_table[c(1:length(peak)),column] <- peak
            }
            peak_table$peakLenth <- NULL
            peak_table_list[[index]] <- peak_table

        } else {

            #### corner case where all identified peaks are length one
            temp <- data.frame(t(peaks))
            colnames(temp) <- paste("peak", seq_along(peaks))
            peak_table_list[[index]] <- temp
            rm(temp)

        }


    }

    names(peak_table_list) <- sample_names
    rm(index)
    return(peak_table_list)
}

