#' @title peakwidth_est
#'
#' @description This function is designed to generate peak width estimates for
#' each TIC peak detected by sliding window analysis.
#'
#' @param peak_vector - numeric vector with names of specific time points of the
#' chromatography data measured. The numeric values correspond to indiciess
#' within the total chromatographic data that span the peak width.
#' @param time - This vector contains the time measurements during the
#' chromatography. This vector is used to match the values in peak_vector
#' to the names in the intensity vector.
#' @param intensity - measured intensity values for chromatorgraphy
#' @param start - numeric index indicating where peak starts. Leave null.
#' @param end - same as above, leave null.
#' @param old_r2 - previous fit of model used to judge recursion of fit.
#'
#' @details This function takes in one peak vector at a time and runs a linear
#' model on the selected start and end points of a peak. By measuring the
#' change of the fit of the model, the function returns an index of values
#' corresponding to a peak. This function works recursively to estimate
#' the width of the peak. Ultimately, it returns the names of the final
#' points in the peak.
#'
#' @return This function returns a scalar value representing the estimated peak
#' width for a given peak.
#'
#' @export
peakwidth_est <- function(peak_vector,
                          time,
                          intensity,
                          start = NULL,
                          end = NULL,
                          old_r2 = NULL) {


  killSwitch <- F

  # check to make sure input values come from vector
  if(!is.numeric(peak_vector)) {
    warning("A non numeric vector was given to peakwidth_est(). This is incorrect. Check the function input.")
  }

  # updating data values to put into regression

  if(!is.null(start)) { # case where we are in second + iteration of algorithm

    end <- end + 1

  } else { # case where the algorithm is run for the first time

    peak_index <- which(time %in% peak_vector)
    if(length(peak_index) == 0) {
      stop("The peak entered here could not be matched to the chromatography data.")
    }

    start <- peak_index[1] - 1
    end <- peak_index[length(peak_index)]

  }

  # terms for lm
  points <- c(start,start-1,start-2,start-3,end,end+1,end+2,end+3)
  intensityObs <- intensity[points]
  modelIndex <- 1:length(intensityObs)

  # correcting na formation within intensity observation vector
  if(any(is.na(intensityObs))) {

    naObsIndex <- which(is.na(intensityObs))
    modelIndex <- modelIndex[-naObsIndex]
    intensityObs <- intensityObs[-naObsIndex]
    killSwitch <- T

  }

  # running smoothing spline on the data
  splineOut <- smooth.spline(modelIndex, intensityObs)
  splineObs <- splineOut$fit$coef
  splineIndex <- 1:length(splineObs)

  # running a linear model on the outcome of the spline
  chrom_table <- data.frame(splineIndex, splineObs)
  model <- lm(formula = splineObs ~ splineIndex, data = chrom_table)
  new_r2 <- summary(model)$r.squared

  # used for comparison with previous itterations
  if(is.null(old_r2)) {
    old_r2 <- 0
  }

  # recursive case - returns previously calculated model fit if improvement
  if((old_r2 > new_r2 & old_r2 > .75) | 0.05 >= (1 - old_r2) | killSwitch == T) {

    # make sure to return numerical index of fit - not the values being compared
    peak_width <- c(peakStart = start-3, peakEnd = end+3)
    return(peak_width)

  } else {
    peakwidth_est(peak_vector, time, intensity, start, end, old_r2 = new_r2)
  }
}

