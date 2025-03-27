data <- read.csv("/scratch/dbihnam/lsc585/turtleProject/alignedBam/flagstatReports/flagstat_summary.csv")

data$propPairedPercent <- (data$Properly.Paired.Reads / data$Total.Reads) * 100

dataGC <- read.csv("/scratch/dbihnam/lsc585/turtleProject/alignedBam/GC_test/flagstatReports/GCflagstat_summary.csv", header = FALSE)

dataGC$propPairedPercent <- (dataGC$V3 / dataGC$V2) * 100
