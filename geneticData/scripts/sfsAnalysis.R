#This script is intended to plot the site frequency spectrum data generated
#by running sfs.sh/sfsLocation.sh in regular and log-transformed form

#Load libraries-----------------------------------------------------------------
library(ggplot2)

#Data wrangling-----------------------------------------------------------------
#Load the frequency data
freqData <- read.table("/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/sfs/sfs_preprocessed.frq", header = TRUE, comment.char = "")

#Calculate the number of alternate alleles (NUM_ALT)
freqData$NUM_ALT <- round(freqData$FREQ. * freqData$N_CHR)

#Calculate the minor allele count (folded SFS)
freqData$MINOR_COUNT <- pmin(freqData$NUM_ALT, freqData$N_CHR - freqData$NUM_ALT)

#Create the SFS by counting the number of SNPs for each minor allele count
maxCount <- max(freqData$MINOR_COUNT, na.rm = TRUE)
sfs <- table(factor(freqData$MINOR_COUNT, levels = 0:maxCount))

#Convert to a data frame for plotting
sfsDf <- data.frame(
  Allele_Count = as.numeric(names(sfs)),
  Frequency = as.numeric(sfs)
)

#Remove the k=0 bin (if any monomorphic sites are in the dataset)
sfsDf <- sfsDf[sfsDf$Allele_Count > 0, ]

#Plotting-----------------------------------------------------------------------
#Plot regular SFS
ggplot(sfsDf, aes(x = Allele_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#01796F", color = "white", size = 0.2) +
  theme_minimal() +
  labs(x = "Minor Allele Count (k)", y = "Number of SNPs", 
       title = expression("Site Frequency Spectrum (Folded) for " * italic("Chelonia mydas ") * "(n=237)")) +
  scale_x_continuous(breaks = seq(0, max(sfsDf$Allele_Count), by = 20)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#Log-transformed SFS
ggplot(sfsDf, aes(x = Allele_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#01796F", color = "white", size = 0.2) +
  theme_minimal() +
  labs(x = "Minor Allele Count (k)", y = "Number of SNPs (log scale)", 
       title = expression("Log-Transformed Site Frequency Spectrum (Folded) for " * italic("Chelonia mydas ") * "(n=237)")) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0, max(sfsDf$Allele_Count), by = 20)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
