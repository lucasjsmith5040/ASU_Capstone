# Read the PCA eigenvector file
pca <- read.table("/scratch/dbihnam/lsc585/turtleProject/variants/filtered/merged237/merged_pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:20))

# Optionally, merge with your metadata
metadata <- read.csv("sample_metadata.csv")  # Your file with group info, e.g., age or location
pca_data <- merge(pca, metadata, by="IID")

# Plot the first two PCs
library(ggplot2)
ggplot(pca_data, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=2) +
  theme_minimal() +
  labs(title="PCA of Merged VCF Data",
       x="Principal Component 1",
       y="Principal Component 2")
