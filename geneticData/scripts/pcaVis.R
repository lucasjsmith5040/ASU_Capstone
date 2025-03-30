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

#Set color palette
library(RColorBrewer)
colorPalette <- brewer.pal(min(length(unique(mergedData$BinnedLocation)), 12), "Paired")

#Plot instructions
fig <- plot_ly(data = mergedData, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~BinnedLocation,  # Color based on BinnedLocation
               colors = colorPalette,   # Apply custom color palette
               type = 'scatter3d', 
               mode = 'markers',
               size = 2.5)  # Adjust marker size

#Plot axis/labels
fig <- fig %>% layout(
  title = paste("PCA: Genomic SNPs Identified in Eastern Mediterranean <i>Chelonia mydas</i> samples (n=237)", sep=""),
  scene = list(
    xaxis = list(title = 'PC1'),
    yaxis = list(title = 'PC2'),
    zaxis = list(title = 'PC3')
  )
)

# Show the plot
fig








