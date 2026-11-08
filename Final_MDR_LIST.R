# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(rlang)
install.packages("writexl")
library(writexl)
library(openxlsx)


## setwd()
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/")

library(readxl)

# Assuming Correlation.xlsx is in the current working directory and you want to read the first sheet
data <- read_excel("Farm_Calf_Ecoli.xlsx", sheet = 7)

# Define the cut-off values for each antibiotic (replace with your actual cut-off values)
cut_offs <- c(
  AUG = 13.01,
  CAZ = 17.01,
  CHL = 12.01,
  CIP = 21.01,
  CTX = 22.01,
  STR = 11.01,
  TET = 11.01,
  TS = 10.01
)

# Create a new column to count the number of resistances
data <- data %>%
  rowwise() %>%
  mutate(ResistanceCount = sum(c_across(names(cut_offs)) < unlist(cut_offs)))

# Filter isolates that are resistant to 3 or more antibiotics
resistant_3_or_more <- data %>%
  filter(ResistanceCount >= 3)

# Select the columns for output (Isolate, Farm, and antibiotic resistance zones)
final_output <- resistant_3_or_more %>%
  select(IsolateID, FarmID, everything())

# Print the result
print(final_output)

View(final_output)

#Save the result into a new Excel file
output_file_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/resistant_isolates_highlighted_25.xlsx"
write_xlsx(highlighted_resistance, output_file_path)
