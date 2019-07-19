#' @title estimate_loess
#'
#' @param xs - an xcmsSet object following peakpicking and grouping
#' @param sampleClasses - String vector representing the classes each sample
#' belogns to.
#' @param max_missing - Estimate for maximum missing peaks within a group to
#' be considered for retention time correction.
#' @param threshold - peak width parameter to estimate retention time
#' correction does not overfit the data.
#'
#' @return This function will return the maximum span that will still
#' provide retention time correction across features over time scales
#' as large as the estimated peak dispersion parameters.
#'
#' @export
estimate_loess <- function(xs, max_missing, sampleClasses, threshold) {

    # reactive variable input that could be caryover from the previous reactive
    ##input
    ## ENTER REACTIVE INPUT HERE TO
    sample_pop <- which(sampleClasses != "blank")

    sample_max_dists <- list()
    itterator <- seq(0.3,0.7,by = 0.01)

    sample_max_dists <- lapply(1:length(itterator), function(i) {
    ii <- itterator[i]
    loess_params <- list(object = xs,
             smooth = "loess",
             missing = max_missing,
             extra = 0,
             span = ii, # itterate with span by visualing data
             family = "gaussian",
             col = xs@phenoData[,"class"])

    suppressWarnings(rc.obi <- do.call(xcms::retcor.peakgroups, loess_params))

    distance <- list()
    for(j in 1:length(rc.obi@rt$raw)) {
        distance[[j]] <- abs(rc.obi@rt$raw[[j]] - rc.obi@rt$correcte[[j]])

    }
    sample_max_dists <- distance[sample_pop] %>%
        lapply(max) %>%
        unlist()
    return(sample_max_dists)

    })

    max_span <- T
    k <- 1
    while(max_span == T && k < length(itterator)) {
    max_span <- any(sample_max_dists[[k]] > threshold)
    k <- k + 1
    }
    estimated_param <- itterator[k]

    return(estimated_param)
}
