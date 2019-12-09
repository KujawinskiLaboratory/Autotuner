context("Performing dataset wide parameter selection")

# The eicParamEsts object was created by running the code within
# inst/makeTestData/
Autotuner <- readRDS(system.file("extdata/AutotunerObj.rds", package="Autotuner"))
eicTable <- readRDS(system.file('extdata/eicParamsEsts.rds',
                                package = "Autotuner"))

returnParams <- returnParams(eicTable, Autotuner)

test_that(desc = "Testing Max Peakwidth Estimate",
          code = {
              expect_equal(class(returnParams), "list")
              expect_equal(length(returnParams), 2)
              expect_equal(dim(returnParams$ticParams), c(3,2))
              expect_equal(dim(returnParams$eicParams), c(7,4))
          })
