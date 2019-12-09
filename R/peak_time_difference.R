#' @title peak_time_difference
#'
#' @description This function is designed to return a data.frame containing info
#' on how
#'
#' @param Autotuner An Autotuner object containiing a table of peak width
#' values extracted with the function peak_width_table.
#'
#' @details This function is designed to determine what are the retention time
#' differences between peaks that are effectively the same between samples.
#' The similarity in peaks is determined by a threshold in retention time
#' similarity between peaks. This function returns the max peak width between
#' samples, and the time difference between peaks across samples in a data frame
#' object. The current and next row indexes are given to go back to the
#' peaktable object to plot peaks.
#'
#' @return This function returns a data.frame of peaks matched over time.
peak_time_difference <- function(Autotuner) {


    peak_table <- getAutoPeak_table(Autotuner)

    # Checking Input ----------------------------------------------------------
    assertthat::assert_that(nrow(peak_table) > 0, msg =
                                paste("No peakwidth table observed.",
                                      "Check output of peakwidth_table."))


    # Sorting peaks across samples by mid point time --------------------------
    peak_table$index <- seq_len(nrow(peak_table))
    sorted_peakwidths <- peak_table[order(peak_table$Mid_point_time,
                                          decreasing = TRUE),]

    ## makes the assumption that dispersion is not as bad as peak distance
    ## from one another
    ## makes sense - if peaks are very far from one another
    notZero <- function(numVector) {
        return(numVector[numVector != 0])
    }


    # Generating a pseudo rt correction threshold  ----------------------------
    ## this is the same as in the "extract_peaks" function
    ## using median as the soft threshold since it is statistically robust.
    ## Now estimating it from all ordered quantities within the data.
    threshold <- peak_table
    threshold <- split(x = threshold, f = peak_table$Sample)
    threshold <- lapply(X = threshold, FUN = function(df) {

        if(nrow(df) > 1) {
            startDiff <- sort(notZero(df$Start_time))
            startDiff <- diff(startDiff)
            startDiff <- min(startDiff)
            endDiff <- sort(notZero(df$End_time))
            endDiff <- diff(endDiff)
            endDiff <- min(endDiff)
            maxDiff <- sort(notZero(df$Maxima_time))
            maxDiff <- diff(maxDiff)
            maxDiff <- min(maxDiff)
        } else {
            startDiff <- notZero(df$Start_time)
            endDiff <- notZero(df$End_time)
            maxDiff <- notZero(df$Maxima_time)
        }

        return(list(startDiff, endDiff, maxDiff))
    })
    threshold <- unlist(threshold)
    threshold <- median(threshold)

    # initializing storage objects --------------------------------------------
    ## initializing storage data
    current_row <- 1
    next_row <- current_row + 1
    storage_counter <- 1
    matchingPeaks <- data.frame()



    # Matching Peaks between Samples ------------------------------------------
    while(next_row <= nrow(sorted_peakwidths)) {


        # checking sample origin -----------------------------------------------
        ## Case 1 - adj peaks come from the same sample
        if(identical(sorted_peakwidths$Sample[current_row],
                     sorted_peakwidths$Sample[next_row])) {

            ## make sure the next sample comes from a different data samples
            current_row <- next_row
            next_row <- current_row + 1

        } else {

        ## case 2 - they come from different samples

        # Checking peak identity -----------------------------------------------
        ## checking if the peak is the same across samples
        ## condition to elimitate two overlapping peaks
        # either the midpoint or the maxima of both peaks have to be within the
        # interval of the minimum start and maximum end points of the data

            row_subset <- c(current_row,next_row)

            interval <- c(start = max(sorted_peakwidths$Start_time[row_subset]),
            end = min(sorted_peakwidths$End_time[row_subset]))

            peakDiffStart <- abs(diff(sorted_peakwidths$Start_time[row_subset]))
            peakDiffEnd <- abs(diff(sorted_peakwidths$End_time[row_subset]))
            peakDiffMid <- abs(diff(
                sorted_peakwidths$Mid_point_time[row_subset]))
            peakDiffMax <- abs(diff(sorted_peakwidths$Maxima_time[row_subset]))

            intervalVals <- c(sorted_peakwidths$Mid_point_time[row_subset],
            sorted_peakwidths$Maxima_time[row_subset])

            ## makes sure peaks are not losely overlapping - increase thresohld
            intervalCheck <- all(interval["start"] < intervalVals &
                                     interval["end"] > intervalVals)
            thresholdCheck <- peakDiffStart <= threshold &
                peakDiffEnd <= threshold

            ## Case 1 - peaks match from the set threshold
            if(all(thresholdCheck,intervalCheck)) {

            storeData <- data.frame(
                Max_width = max(sorted_peakwidths$peak_width[row_subset],
                                na.rm = TRUE),
                Start_diff = peakDiffStart,
                End_diff = peakDiffEnd,
                Mid_diff = peakDiffMid,
                Max_diff = peakDiffMax,
                cur_row = sorted_peakwidths$index[current_row],
                next_row = sorted_peakwidths$index[next_row],
                max_intensity = max(sorted_peakwidths$Max_intensity[row_subset])
            )

            matchingPeaks <- rbind(matchingPeaks, storeData)

            # only updating the next value
            next_row <- next_row + 1
            storage_counter <- storage_counter + 1

            ## Case 2 - peaks are different from one another.
        } else {

            current_row <- next_row
            next_row <- current_row + 1
        }

        }
    } # End of while loop


    # Checking performance of peak matching -----------------------------------
    if(nrow(matchingPeaks) == 0) {
        stop(paste("No peaks were matched across samples.",
                   "Consider asking extract_peaks to return more peaks."))
    }
    matchingPeaks$index <- 0


    # initializing second step input variables --------------------------------
    current_row <- 1
    next_row <- current_row + 1
    # Identifying common peaks in 2 or more rows ------------------------------
    while(next_row <= nrow(matchingPeaks)) {

        ## case 1 - Grouping matches between more than 2 peaks
        if(matchingPeaks$cur_row[current_row] ==
           matchingPeaks$cur_row[next_row]) {

            rowIndex <- c(current_row,next_row)
            matchingPeaks$max_intensity[rowIndex] <-
            max(matchingPeaks$max_intensity[rowIndex])

            current_row <- next_row + 1
            next_row <- current_row + 1

        } else {

            current_row <- next_row
            next_row <- next_row + 1

        }
    }


    # organizing and exporting matches ----------------------------------------
    index_vector <- matchingPeaks$cur_row
    index_vector <- unlist(index_vector)
    index_vector <- rle(index_vector)

    temp <- lapply(seq_along(index_vector$lengths), function(i) {
        rep(i, index_vector$lengths[i])
    })
    matchingPeaks$index <- unlist(temp)
    rm(temp)

    return(matchingPeaks)
}

