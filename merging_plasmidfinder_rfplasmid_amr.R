# merge_plasmid_prediction.R
# Joins AMRFinder/PlasmidFinder results with plasmid prediction results
# Matches: file1$SEQUENCE <-> file2$contigID
# Output columns: #FILE, prediction, votes chromosomal, votes plasmid, contigID, GENE

library(dplyr)
library(readr)
library(readxl)

# ── File paths ────────────────────────────────────────────────────────────────
# Update these paths to point to your actual files
file1_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/PLASMIDFINDER.xlsx"   # PlasmidFinder / AMRFinder results (TSV)
file2_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/prediction_Ecoli_plasmids.csv"   # Plasmid prediction results (TSV)
output_path <- "merged_Ecoli_plasmids.tsv"

# ── Read files ────────────────────────────────────────────────────────────────
# File 1: tab-separated, header starts with '#FILE'
df1 <- read_excel(
  file1_path)

# Rename '#FILE' to 'FILE' so dplyr handles it cleanly
names(df1)[names(df1) == "#FILE"] <- "FILE"

# File 2: tab-separated
df2 <- read_csv(
  file2_path,
  col_types = cols(.default = "c")
)

# ── Join on SEQUENCE == contigID ──────────────────────────────────────────────
merged <- inner_join(
  df1 %>% select(FILE, SEQUENCE, GENE),
  df2 %>% select(prediction, `votes chromosomal`, `votes plasmid`, contigID),
  by = c("SEQUENCE" = "contigID")
) %>%
  select(
    `#FILE`              = FILE,
    prediction,
    `votes chromosomal`,
    `votes plasmid`,
    contigID             = SEQUENCE,
    GENE
  )

# ── Write output ──────────────────────────────────────────────────────────────
write_tsv(merged, output_path)

cat(sprintf(
  "Done. %d matching rows written to '%s'.\n",
  nrow(merged), output_path
))

setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/")

file3_path   <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/amrfinder_combined_2.xlsx"  # <- update filename as needed
step2_output <- "merged_Ecoli_plasmids_AMR.tsv"

# If running Step 2 separately, load the Step 1 output instead:
merged1 <- read_tsv(output_path, col_types = cols(.default = "c"))

df3 <- read_excel(file3_path, col_types = "text")
names(df3) <- trimws(names(df3))

merged2 <- inner_join(
  merged1,
  df3 %>% select(`Contig id`, `Gene symbol`, `Sequence name`, `Element subtype`, Class),
  by = c("contigID" = "Contig id")
) %>%
  select(
    `#FILE`,
    prediction,
    `votes chromosomal`,
    `votes plasmid`,
    contigID,
    GENE,
    `Gene symbol`,
    `Sequence name`,
    `Element subtype`,
    Class
  )

write_tsv(merged2, step2_output)
cat(sprintf("Step 2 done. %d rows written to '%s'.\n", nrow(merged2), step2_output))
