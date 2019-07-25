
# itterating between samples ----------------------------------------------
Autotuner <- readRDS(system.file("extdata/Autotuner.rds", package="Autotuner"))


currentTable <- peak_table[Autotuner@peak_table$Sample == 1,]
currentFile <- Autotuner@file_paths[1]

    # Adding msnbase functionality to replace mzR API
msnObj <- suppressMessages(MSnbase::readMSData(files = currentFile,
                                               mode = "onDisk",
                                               msLevel. = 1))

header <- suppressWarnings( MSnbase::header(msnObj))

#saveRDS(header, file = here::here("inst/extdata/header.rds"))

allMzs <- MSnbase::mz(msnObj)
allInt <- MSnbase::intensity(msnObj)

mzDb <- list()
for(i in seq_along(allInt)) {
    mzDb[[i]] <- cbind(mz = allMzs[[i]],
                       intensity = allInt[[i]])
}
rm(allMzs, allInt, msnObj, i)

    # going through each peak from a sample -----------------------------
curPeak <- 1

start <- currentTable[curPeak,"Start_time"]
end <- currentTable[curPeak,"End_time"]
width <- currentTable$peak_width[curPeak]
observedPeak <- list(start = start, end = end)

usethis::use_data(mzDb, header, observedPeak)

