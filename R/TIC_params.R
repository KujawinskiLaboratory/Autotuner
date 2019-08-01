#' @title TIC_params
#'
#' @description This function is designed to return the parameters related to
#' chromatography gathered during the TIC peak selection steps. An estimate is
#' given for maximum and minimum peak width as well as the bandwidth parameter
#' used in grouping.
#'
#' @param peak_difference A data.frame containing information on retention
#' time differences between peaks.
#' @param peak_table A data.frame containing information on peak width values
#' extracted with the function peak_width_table.
#'
#' @return Returns a set of parameters to run xcms.
TIC_params <- function(peak_table, peak_difference) {

    # min width estimation
    min_width <- min(peak_table$peak_width)

    # max width estimation
    # Consider ways of increasing this wrt chrom types
    max_width <- max(peak_table$peak_width, na.rm = TRUE)

    group_diff <- max(peak_difference$Mid_diff, na.rm = TRUE)

    return(list(max_width = max_width,
                min_width = min_width,
                group_diff = group_diff))
}
