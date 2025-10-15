# Peak processing of recalibrated datasets

SHARING <- FALSE

library("Cardinal")
RNGkind("L'Ecuyer-CMRG")

setCardinalBPPARAM(MulticoreParam(progressbar=TRUE))
setCardinalVerbose(TRUE)


# Directories and filenames
DATAPATH <- "/Users/lakkimsetty.s/Documents/9-data"
CURRPATH <- "14-Buck-Institute/imzMLs/processed-in-R/recalibrated-images"
COND <- "CAD" ####
IND <- "6919" ####
ID <- "b" ####
SLICE <- "lateral" ####

PREFIX <- paste0(IND, ID, "-", COND, "-", SLICE, "-")
mse_filename <- paste0(PREFIX, "calibrated.imzML")


# mse
mse_calibrated <- readMSIData(file.path(DATAPATH, CURRPATH, COND,
    IND, SLICE, mse_filename))

# Peak picking was done similar to the above for all images and 
# a common peaklist was created. Then I repicked the peaks with the 
# common peaklist as ref below. The common peaklist is also included 
# the archive. 

mse_pp <- mse_calibrated |>
    peakPick(method="diff", SNR=7) |>
    peakAlign(ref=mz_pp, tolerance=80, units="ppm") |> 
    peakFilter(freq.min=0.01) |> 
    process() 

writeMSIData(mse_pp, file.path(DATAPATH,
                               "14-Buck-Institute/imzMLs/processed-in-R",
                               "peakpicked-images", COND, IND, SLICE,
                               paste0(PREFIX, "peakpicked.imzML")))

gc()
all.equal(mz_pp, mz(mse_pp))


