install.packages("readxl")
install.packages("tidyverse")
install.packages("reshape2")

# Load packages
library(readxl)
library(tidyverse)
library(reshape2)
library(pheatmap)
## setwd()
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/")

# Assuming Correlation.xlsx is in the current working directory and you want to read the first sheet
df <- read_excel("MDR_list_final.xlsx", sheet = 4)

# Define resistance cut-off values
cut_offs <- c(
  AUG = 13.01,
  CAZ = 17.01,
  CHL = 12.01,
  CIP = 21.01,
  CTX = 22.01,
  STR = 11.01,
  TET = 11.01,
  TS  = 10.01
)

# Extract the relevant columns (IsolateID + antibiotics)
ab_cols <- names(cut_offs)
df_ab <- df %>% select(IsolateID, all_of(ab_cols))

# Convert to long format
long_df <- pivot_longer(df_ab, cols = -IsolateID, names_to = "Antibiotic", values_to = "Zone")

# Classify as Resistant or Sensitive
long_df <- long_df %>%
  mutate(
    Resistance = ifelse(Zone < cut_offs[Antibiotic], "R", "S")
  )

# Plot the heatmap
ggplot(long_df, aes(x = Antibiotic, y = IsolateID, fill = Resistance)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("S" = "#91d18b", "R" = "#f05d5e")) +
  labs(
    title = "Antibiotic Resistance Heatmap",
    x = "Antibiotic",
    y = "Isolate ID"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))  # for tight labels_

##NOW WITH A CLUSTERING OF THE Y AXIS

# Convert to a resistance matrix: R = 1, S = 0
res_matrix <- df_ab
for (ab in ab_cols) {
  res_matrix[[ab]] <- ifelse(res_matrix[[ab]] < cut_offs[ab], 1, 0)
}
row.names(res_matrix) <- res_matrix$IsolateID
res_matrix <- res_matrix %>% select(-IsolateID)

# Plot heatmap with clustering on Y-axis (rows)
pheatmap(
  mat = res_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # you can set to TRUE if you want column clustering too
  color = c("green", "red"),
  legend_breaks = c(0, 1),
  legend_labels = c("Sensitive", "Resistant"),
  main = "Resistance Heatmap (Clustered by Isolate)",
  fontsize_row = 6,
  fontsize_col = 10,
  border_color = "grey90"
)


