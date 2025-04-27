##This R script is intended to plot a mutational signature plot mimicing the style
##used in COSMIC mutational signatures

##The input file for this analysis is generated from the python script:
##extractMutSigs.py

#Load libraries-----------------------------------------------------------------
library(dplyr)
library(ggplot2)

#Data wrangling-----------------------------------------------------------------
#Load in data
data <- read.csv("/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/mutationalSignatures/captiveSignatures.csv")  # Replace with your CSV file path

#Extract mutation type and group by mutationType
data$mutationType <- sapply(data$Mutation, function(x) {
  bases <- unlist(strsplit(x, ">"))
  ref_base <- substr(bases[1], 2, 2)
  alt_base <- substr(bases[2], 2, 2)
  paste0(ref_base, ">", alt_base)})

#Create values of all possible mutation types
bases <- c("A", "C", "G", "T")
sbsType <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

#Create all possible trinucleotide contexts for C and T
allContext <- expand.grid(
  fivePrime = bases,
  ref = c("C", "T"),
  threePrime = bases)

#Filter and create mutation strings for each substitution type
allMutation <- data.frame()
for (sub in sbsType) {
  #C or T
  ref_base <- substr(sub, 1, 1)
  #A, G, T, etc.
  alt_base <- substr(sub, 3, 3)
  contexts <- allContext[allContext$ref == ref_base, ]
  contexts$Mutation <- paste0(contexts$fivePrime, ref_base, contexts$threePrime, ">", contexts$fivePrime, alt_base, contexts$threePrime)
  contexts$mutationType <- sub
  allMutation <- rbind(allMutation, contexts)}

#Merge data with all 96 contexts, fill in any zeros
allMutation <- allMutation[, c("Mutation", "mutationType")]
mergedData <- left_join(allMutation, data[, c("Mutation", "Count", "mutationType")], by = c("Mutation", "mutationType"))
mergedData$Count[is.na(mergedData$Count)] <- 0

#Convert to percentages
totalCount <- sum(data$Count)
mergedData$Percentage <- (mergedData$Count / totalCount) * 100

#Plot---------------------------------------------------------------------------
ggplot(mergedData, aes(x = Mutation, y = Percentage, fill = mutationType)) +
  geom_bar(stat = "identity", color = "white", width = 0.7) +
  facet_wrap(~ mutationType, scales = "free_x", nrow = 1) +
  labs(title = "Captive: Mutational Signature",
       y = "Mutation Frequency (%)",
       x = "Trinucleotide Mutational Context") +
  scale_y_continuous(limits = c(0, 5)) +  # Adjusted for percent
  scale_fill_manual(values = c("C>A" = "#1dbdf1", "C>G" = "#060709", "C>T" = "#e72626",
                               "T>A" = "#cacaca", "T>C" = "#a0cf65", "T>G" = "#edc7c4")) +
  theme_minimal(base_family = "Calibri") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank())
