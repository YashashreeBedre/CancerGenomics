setwd("C:/Users/YashashreeBedre/OneDrive/Desktop/cancer genomics")

data<- read.table("GSE200146_counts.Cell.Lines.txt", 
                  header=TRUE, row.names=1, sep = '\t')
data= data[,-1]
head(data)


sizefactors= colSums(data)
columnsums = apply(data, 2, sum)


normalizedMatrix <- apply(data, 2, function(x)(x/sum(x))*1000000)
normalizedMatrixlog= log2(normalizedMatrix+1)
dim(normalizedMatrixlog)

# Assuming you have a numeric matrix named "data"
standardizedData <- matrix(NA, nrow = nrow(data), ncol = ncol(data)) 
# Create an empty matrix for standardized values
dim(standardizedData)

for (i in 1:100) {
  for (j in 1:100) {
    mean_value <- mean(normalizedMatrixlog[i, ])
    sd_value <- sd(normalizedMatrixlog[i, ])
    standardizedData[i, j] <- (normalizedMatrixlog[i, j] - mean_value) / sd_value
  }
}
# Calculate Counts Per Million (CPM)
cpm <- apply(data, 2, function(x) (x / sum(x) * 10^6))
dim(a)

# Define a log transformation function
logtransform <- function(cpm) {
  cpm1 <- log2(cpm + 1)
  return(cpm1)
}
# Apply log transformation to the CPM data
cpm2 <- logtransform(cpm)
dummy <- cpm2


mean_value <- mean(normalizedMatrixlog)
sd_value <- sd(normalizedMatrixlog)
z_scores <- (normalizedMatrixlog - mean_value) / sd_value

library(ComplexHeatmap)
library(circlize)
Heatmap(z_scores[1:100, ], col=colorRamp2(c(-2, 0, 2),c('orange', 'white', 'purple')))
saveRDS(z_scores, "z_scores.rds")

# Calculate variance for each gene
gene_var <- apply(data[1:100], 1, var)
print(gene_var)

# Task 2: Sort the total variance in decreasing order and select the top 100 rows
new_variance <- sort(gene_var, decreasing = TRUE)
sel <- names(new_variance)[1:100]
cpm_new <- cpm2[sel,]

# Task 3: Calculate the Z scores for the selected rows
calculateZ_SCORE1 <- function(cpm_new) {
  Tmean1 <- apply(cpm_new, 1, mean)
  tsd1 <- apply(cpm_new, 1, sd)
  z_score1 <- cpm_new
  for (i in 1:nrow(cpm_new)) {
    z_score1[i, 1:ncol(z_score1)] <- (cpm_new[i, 1:ncol(cpm_new)] - Tmean1[i]) / tsd1[i]
  }
  return(z_score1)
}

zscore_new <- calculateZ_SCORE1(cpm_new)

# Task 4: Create the heatmap for the Z scores
Heatmap(zscore_new, col = colorRamp2(c(-2, 0, 2), colors = c("orange", "white", "purple")))

# Load annotation data
anno <- read.csv("C:/Users/YashashreeBedre/Downloads/meta.csv")
names(anno) 

# Color scales for Age, Gender, and X columns
col1 <- list(call = colorRamp2(c(-2, 0, 2), colors = c("orange", "white", "purple")),
             agee = colorRamp2(c(-2, 0, 2), colors = c("orange", "white", "purple")),
             Gender = c("Female" = "red", "Male" = "blue"))

# Annotation data
ann <- data.frame(call = anno$x, Age = anno$Age, Gender = anno$Gender)

# Create heatmap annotation for Age, Gender, and X with color scales
ha <- HeatmapAnnotation(df = ann, col = col1)

# Create the heatmap with annotations and legends
pdf(file= "C:/Users/YashashreeBedre/OneDrive/Desktop/cancer genomics/myplot.pdf", width = 35, height = 15)
Heatmap(zscore_new, top_annotation = ha)

# Save the heatmap with annotations and legends
dev.off()
