## This R script is intended to calculate the read mapping percentage of all BAM
## flagstat reports, and visualize it in a dataframe

## This was mainly useful in determining a cutoff % that would maintain study size
## AND sample quality in downstream analysis

# Read in parsed flagstat document with merged flagstat reports of all samples
data <- read.csv("/geneticData/alignedBam/flagstatReports/flagstat_summary.csv")

# Calculate read mapping percentage
data$propPairedPercent <- (data$Properly.Paired.Reads / data$Total.Reads) * 100

# Read in parsed flagstat document with merged flagstat reports of all samples (with high GC alignment, not used)
dataGC <- read.csv("/geneticData/alignedBam/GC_test/flagstatReports/GCflagstat_summary.csv", header = FALSE)

# Calculate read mapping percentage for high GC alignment (not used)
dataGC$propPairedPercent <- (dataGC$V3 / dataGC$V2) * 100
