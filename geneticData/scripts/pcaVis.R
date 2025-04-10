###NEED TO UPDATE 3D PCA PLOT BEFORE UPDATING GITHUB!!! MAKE SURE TO REPLACE PLOTS IN GITHUB AS WELL!!
#Import libraries
library(ggplot2)
library(dplyr)

#Load pca
pcaData <- read.table("/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/pca/10/pcaResults.eigenvec", header = FALSE, stringsAsFactors = FALSE)
colnames(pcaData) <- c("FID", "IID", paste0("PC", 1:10))

#Extract sample ID and remove file extension
pcaData$IID <- sub(".*/", "", pcaData$IID)
pcaData$IID <- gsub("\\.sorted\\.bam$", "", pcaData$IID)

#Load metadata
metadata <- read.csv("/scratch/dbihnam/lsc585/turtleProject/sequencingMetadata.csv")

#Merge PCA results with metadata
mergedData <- merge(pcaData, metadata, by.x = "IID", by.y = "SequencingSampleID")

#Plot PCA in 2D by sample location
ggplot(mergedData, aes(x = PC1, y = PC2, color = BinnedLocation)) + 
  geom_point(size = 2) +
  labs(title = expression("PCA: Genomic SNPs Identified in Eastern Mediterranean " * italic("Chelonia mydas ") * "samples (n=237)"), x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_discrete(name = "Sample Location")

#Plot PCA in 3D by sample location
#Load the plotly package
#install.packages("plotly")
library(plotly)

#Define the color palette
colors <- c(
  "forestgreen",  # Ashkelon
  "red",          # Captive Breeding
  "darkgray",     # Caesarea
  "blue",         # Central-Northern Coast (Other)
  "cyan",         # Gdor Nature Reserve
  "purple",       # Habonim Beach Nature Reserve
  "orange",       # Hof Galim Megdim
  "chocolate",    # Nahariyya-Betzet
  "darkred",      # Netanya
  "hotpink",      # Nitzanim
  "gold",         # Stranded
  "magenta",      # Yavne
  "limegreen"     # Zikim
)

#Ensure the colors match the order of BinnedLocation levels
unique_locations <- levels(factor(mergedData$BinnedLocation))
if (length(unique_locations) != length(colors)) {
  stop("Number of unique BinnedLocation levels does not match the number of colors!")
}
names(colors) <- unique_locations

#Plot instructions
fig <- plot_ly(data = mergedData, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~BinnedLocation,
               colors = colors,
               type = 'scatter3d', 
               mode = 'markers',
               marker = list(size = 3))

#Plot axis/labels and adjust legend
fig <- fig %>% layout(
  title = paste("PCA: Genomic SNPs Identified in Eastern Mediterranean <i>Chelonia mydas</i> samples \n(n=237, 44.33% of Variance Explained)", sep=""),
  titlefont = list(size = 18),
  scene = list(
    xaxis = list(title = 'PC1 - 18.45% Variance'),
    yaxis = list(title = 'PC2 - 14.09% Variance'),
    zaxis = list(title = 'PC3 - 11.79% Variance')
  ),
  legend = list(
    title = list(text = "Sample Location", font = list(size = 18)),
    font = list(size = 16),
    itemsizing = "constant"
  )
)

#View plot
fig





