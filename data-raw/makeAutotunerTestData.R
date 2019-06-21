
# itterating between samples ----------------------------------------------
system.file("peak_table.rda", package = "Autotuner")
data("peak_table", package="Autotuner")
data("peak_difference", package="Autotuner")
data("Autotuner", package="Autotuner")

currentTable <- peak_table[peak_table$Sample == 1,]
currentFile <- Autotuner@file_paths[1]

    # Adding msnbase functionality to replace mzR API
msnObj <- suppressMessages(MSnbase::readMSData(files = currentFile,
                                               mode = "onDisk",
                                               msLevel. = 1))

header <- suppressWarnings( MSnbase::header(msnObj))
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

