library(readxl)
library(tidyverse)
library(vegan)      # for vegdist (Jaccard)
library(ggdendro)   # for ggplot-style dendrograms
# or library(ggtree) if you want to link it to a metadata strip like your other trees
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/")
# 1. Read in the data
vir <- read_excel("Virulence_Ecoli.xlsx", sheet = 1)

# 1. Build presence/absence matrix: rows = samples, columns = virulence genes
vir_wide <- vir %>%
  distinct(sample, `Virulence factor`) %>%
  mutate(Presence = 1) %>%
  pivot_wider(
    names_from = `Virulence factor`,
    values_from = Presence,
    values_fill = 0
  )
# Convert to matrix with sample names as rownames
vir_mat <- vir_wide %>%
  column_to_rownames("sample") %>%
  as.matrix()

# 3. Compute Jaccard distance
# vegan::vegdist with method = "jaccard" and binary = TRUE treats data as presence/absence
jaccard_dist <- vegdist(vir_mat, method = "jaccard", binary = TRUE)

# 4. Hierarchical clustering
hc <- hclust(jaccard_dist, method = "average")  # UPGMA is standard for Jaccard-based gene trees

# 5. Plot as a simple dendrogram
plot(hc, hang = -1, cex = 0.7, main = "Hierarchical clustering (Jaccard distance)")

##Now let's export it: 
library(ape)      # for write.tree()

# 1. Export the Newick tree for iTOL
# hclust needs to be converted to a "phylo" object first
tree <- as.phylo(hc)
write.tree(tree, file = "virulence_jaccard_tree.nwk")

# 2. Export the presence/absence matrix as CSV
vir_wide %>%
  write_csv("virulence_presence_absence.csv")

# Optional: also save as xlsx if you want it Excel-formatted
library(writexl)
write_xlsx(vir_wide, "virulence_presence_absence.xlsx")

