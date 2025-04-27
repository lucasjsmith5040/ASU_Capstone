##This R script is intended to take all mutational signatures per region
##and run non-negative matrix factorization on them in order to decompose
##overall signatures into smaller components

#Load in libraries--------------------------------------------------------------
library(dplyr)
library(readr)
library(purrr)
library(NMF)
library(tidyverse)

#Data wrangling-----------------------------------------------------------------
#Navigate to working directory so I don't have to type full paths over and over
setwd("/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/mutationalSignatures/")

#.csv file names
regionFiles <- list(
  Ashkelon = "ashkelonSignatures.csv",
  Captive = "captiveSignatures.csv",
  Ceasarea = "ceasareaSignatures.csv",
  Gdor_Nature_Reserve = "gdorSignatures.csv",
  Habonim_Beach_Nature_Reserve = "habonimSignatures.csv",
  Hof_Galim_Megadim = "hofSignatures.csv",
  Nahariyya_Betzet = "nahariyyaSignatures.csv",
  Netanya = "netanyaSignatures.csv",
  Nitzanim = "nitzanimSignatures.csv",
  Central_Northern_Coast_Other = "otherSignatures.csv",
  Stranded = "strandedSignatures.csv",
  Yavne = "yavneSignatures.csv",
  Zikim = "zikimSignatures.csv"
)

#Define mutation types and trinucleotide contexts
bases <- c("A", "C", "G", "T")
sbsType <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
allContext <- expand.grid(fivePrime = bases, ref = c("C", "T"), threePrime = bases)

#Assemble all 96 trinucleotide mutations
allMutation <- data.frame()
for (sub in sbsType) {
  ref_base <- substr(sub, 1, 1)
  alt_base <- substr(sub, 3, 3)
  contexts <- allContext[allContext$ref == ref_base, ]
  contexts$Mutation <- paste0(contexts$fivePrime, ref_base, contexts$threePrime, ">", contexts$fivePrime, alt_base, contexts$threePrime)
  contexts$mutationType <- sub
  allMutation <- rbind(allMutation, contexts)}
allMutation <- allMutation[, c("Mutation", "mutationType")]

#Function to process each region and create all 96 possible contexts
processRegion <- function(regionName, fileName) {
  filePath <- file.path(getwd(), fileName)
  print(paste("Reading file for region:", regionName, "at", filePath))
  data <- read_csv(filePath)
  #Create mutationType
  data$mutationType <- sapply(data$Mutation, function(x) {
    bases <- unlist(strsplit(x, ">"))
    ref_base <- substr(bases[1], 2, 2)
    alt_base <- substr(bases[2], 2, 2)
    paste0(ref_base, ">", alt_base)})
  #Merge with allMutation to make sure all contexts are present
  mergedData <- left_join(allMutation, data[, c("Mutation", "Count", "mutationType")], by = c("Mutation", "mutationType"))
  mergedData$Count[is.na(mergedData$Count)] <- 0
  mergedData$Region <- regionName
  mergedData}

#Apply function to all regions
combinedData <- map_dfr(names(regionFiles), function(regionName) {
  fileName <- regionFiles[[regionName]]
  processRegion(regionName, fileName)})

#Switch to wide format
nmfInput <- combinedData %>%
  select(Region, Mutation, Count) %>%
  pivot_wider(names_from = Mutation, values_from = Count, values_fill = 0)

#Create signature matrix
regionNames <- nmfInput$Region
nmfMatrix <- nmfInput %>% select(-Region) %>% as.matrix()
rownames(nmfMatrix) <- regionNames

#NMF----------------------------------------------------------------------------
#Prepare input matrix and handle any empty data
nmfMatrix <- apply(nmfMatrix, 2, as.numeric)
nmfMatrix <- nmfMatrix[rowSums(nmfMatrix) > 0, ]
nmfMatrix <- nmfMatrix[, colSums(nmfMatrix) > 0]

#Run NMF
#Rank is number of signatures to extract
#nrun is number of iterations to run
nmfResult <- nmf(nmfMatrix, rank = 4, nrun = 1000)

#Extract W (basis) matrix and label regions
basisMatrix <- nmfResult@fit@W
rownames(basisMatrix) <- c("Ashkelon", "Captive", "Ceasarea", "Gdor_Nature_Reserve", 
                           "Habonim_Beach_Nature_Reserve", "Hof_Galim_Megadim", 
                           "Nahariyya_Betzet", "Netanya", "Nitzanim", 
                           "Central_Northern_Coast_Other", "Stranded", "Yavne", "Zikim")

#Switch to long format
basisDF_long <- basisMatrix %>%
  as.data.frame() %>%
  rownames_to_column("Region") %>%
  pivot_longer(-Region, names_to = "Signature", values_to = "Contribution")

#Normalize contributions to proportions
basisDF_long_norm <- basisDF_long %>%
  group_by(Region) %>%
  mutate(Proportion = Contribution / sum(Contribution)) %>%
  ungroup()

#Plotting-----------------------------------------------------------------------
#Plot count-based de novo signatures
ggplot(basisDF_long, aes(x = Region, y = Contribution, fill = Signature)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Regional Contributions to De Novo Signatures",
       x = "Region",
       y = "Mutation Count") +
  theme_minimal(base_family = "Calibri") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot proportional de novo signatures
ggplot(basisDF_long_norm, aes(x = Region, y = Proportion, fill = Signature)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "De Novo Signature Decomposition by Region",
    x = "Region",
    y = "Proportion") +
  theme_minimal(base_family = "Calibri") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(family = "Calibri", size = 14),
    axis.title = element_text(family = "Calibri", size = 13),
    legend.title = element_text(family = "Calibri", size = 13),
    legend.text = element_text(family = "Calibri", size = 12)) +
  scale_fill_manual(values = c("V1" = "#ee1825", #Raphael
                               "V2" = "#fe8121", #Michaelangelo
                               "V3" = "#3b8bcd", #Leonardo (my favorite)
                               "V4" = "#8645a2")) #Donatello


