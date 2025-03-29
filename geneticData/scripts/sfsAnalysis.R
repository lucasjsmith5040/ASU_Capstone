#Load ggplot2
library(ggplot2)

#Load the frequency data
freq_data <- read.table("/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/sfs/sfs_preprocessed.frq", header = TRUE, comment.char = "")

#Calculate the number of alternate alleles (NUM_ALT)
freq_data$NUM_ALT <- round(freq_data$FREQ. * freq_data$N_CHR)

#Calculate the minor allele count (folded SFS)
freq_data$MINOR_COUNT <- pmin(freq_data$NUM_ALT, freq_data$N_CHR - freq_data$NUM_ALT)

#Create the SFS by counting the number of SNPs for each minor allele count
#We use factor() to ensure all possible counts (0 to max) are included, even if they have 0 SNPs
max_count <- max(freq_data$MINOR_COUNT, na.rm = TRUE)
sfs <- table(factor(freq_data$MINOR_COUNT, levels = 0:max_count))

#Convert to a data frame for plotting
sfs_df <- data.frame(
  Allele_Count = as.numeric(names(sfs)),
  Frequency = as.numeric(sfs)
)

#Remove the k=0 bin (monomorphic sites, if any)
sfs_df <- sfs_df[sfs_df$Allele_Count > 0, ]

#Plot the SFS as a bar plot
ggplot(sfs_df, aes(x = Allele_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", size = 0.2) +
  theme_minimal() +
  labs(x = "Minor Allele Count (k)", y = "Number of SNPs", 
       title = "Site Frequency Spectrum (Folded) for Chelonia mydas") +
  scale_x_continuous(breaks = seq(0, max(sfs_df$Allele_Count), by = 5)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#Re-plot the SFS with adjusted x-axis breaks
ggplot(sfs_df, aes(x = Allele_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", size = 0.2) +
  theme_minimal() +
  labs(x = "Minor Allele Count (k)", y = "Number of SNPs (log scale)", 
       title = "Site Frequency Spectrum (Folded) for Chelonia mydas") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0, max(sfs_df$Allele_Count), by = 20)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
