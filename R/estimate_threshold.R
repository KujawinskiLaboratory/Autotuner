#' @title estimate_threshold
#'
#' @param distance_vector - This vector contains the retention time differences
#' between individual peaks identified within a single sample TIC.
#'
#' @description AS OF 2019-07-02, THIS FUNCTION IS NO LONGER PART OF THE
#' AUTOTUNER ALGORITHM SINCE IT IS NO LONGER CONSIDERED TO BE NECESSARY.
#' At a high level, this function is used to estimate appropriate
#' thresholds to required to distinguish between different peaks accross
#' samples, and replicated ones that only vary due to retention time deviation.
#' It assumes that duplicate peaks are a relatively short distance from one
#' another relative to those that distinct.
#'
#' @details To estimate threshold, all the measured distances are clustered
#' using k-means clustering with an itteritatively increasing number of
#' clusters. The itteration is complete once there are as many clusters as
#' peaks or when the median of the smallest distance cluster is larger than the
#' difference between the minimum and maximum distance values within that
#' cluster.
#'
#' @return A scalar value made up of the largest value predicted to come from
#' peaks identical peaks plus 1.
#'
#' @export
estimate_threshold <- function(distance_vector) {

    k <- 2
    threshold <- 0

    while(k < length(distance_vector) && threshold == 0) {

        cluster_stats <- kmeans(dist(distance_vector), k)

        smallest_dist <- distance_vector[which(distance_vector %in% min(distance_vector))]
        if(!is.null(names(distance_vector))) {
            min_cluster <- cluster_stats$cluster[which(names(cluster_stats$cluster) %in%
                                                       names(smallest_dist)[1])]
        }

        index <- names(cluster_stats$cluster[which(cluster_stats$cluster == min_cluster)])
        cluster_elements <- distance_vector[which(names(cluster_stats$cluster)
                                                  %in% index)]
        if(max(cluster_elements) - min(cluster_elements) < median(cluster_elements)) {
            threshold <- max(cluster_elements) + 1
        }
        k <- k + 1
    }
    return(threshold)
}
