library(ggplot2)

# Load the data
tajima_data <- read.table("/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/tajimaD/grok/tajimaD_output.Tajima.D", header = TRUE)

# Plot
ggplot(tajima_data, aes(x = BIN_START, y = TajimaD, color = CHROM)) +
  geom_point(size = .5) +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Tajima's D")









