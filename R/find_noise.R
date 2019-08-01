#' @title find_noise
#'
#' @param time_trace A chromatographic intensity trace
#'
#' @return This function determines the indexes within a chromatographic
#' intensity trace that are not noise. It does this by asking if a given point
#' is greater than the average of that point and the two bounding points on
#' either side of it.
find_noise <- function(time_trace) {

    not_noise <- vector()
    for(i in 1:length(time_trace)) {

        if(i == 1 || i == length(time_trace)) {
            not_noise[i] <- 0
            next
        }

    if(i > 5 && i < length(time_trace) - 5) {
            m <- mean(time_trace[(i-5):i])
        } else {
            m <- 0
    }

    if(time_trace[i] - time_trace[i - 1] > m &&
        time_trace[i] - time_trace[i + 1] > m) {
        not_noise[i] <- 1
    } else {
        not_noise[i] <- 0
    }

    }
    return(not_noise)
}
