library(tidyverse)
library(readxl)
# Load your data — adjust separator as needed
df <- read_xlsx("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/amrfinder_combined_2.xlsx")  # or read_csv("your_file.csv")

# Assume your data has at minimum two columns:
# - "Sample" : sample name
# - "Gene"   : gene or mutation name

# If you have additional columns (e.g. coverage, identity), 
# just adjust the select() call below

matrix_df <- df %>%
  select(Sample, Gene) %>%          # keep only what we need
  distinct() %>%                    # remove duplicates
  mutate(present = 1) %>%           # mark presence
  pivot_wider(
    names_from  = Gene,
    values_from = present,
    values_fill = -1                 # absent = -1
  )

# Result: Sample column + one column per gene/mutation
# Save it
write_tsv(matrix_df, "presence_absence_matrix.tsv")
