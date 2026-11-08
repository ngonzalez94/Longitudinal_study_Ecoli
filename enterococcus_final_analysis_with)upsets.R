#####################################################################
# Enterococcus faecalis/faecium - full pipeline:
#   1. Count isolates/samples BEFORE filtering (by farm, by plate)
#   2. Filter to only samples in calves_per_farm_final.xlsx (= final,
#      faecal-only samples that made it into the analysis)
#   3. Count isolates/samples AFTER filtering, and report what was lost
#   4. Resistance analysis (isolate + sample level) on the filtered set
#   5. Isolate-level UpSet plot of resistance profiles, split by farm
#####################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ComplexUpset)

setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/")

# ------------------------------------------------------------------
# Species colour palette (E. faecalis light blue, E. faecium light pink)
# ------------------------------------------------------------------
species_colors <- c(
  "faecalis" = "#AEE1F0",
  "faecium"  = "#F7C6C7"
)

# Farm colour palette (update names to match your actual FarmID values -
# check with sort(unique(entero$FarmID)) after building `entero` below)
farm_colors <- c(
  "1" = "#EEC240",
  "2" = "#2C77B2",
  "3" = "#C80012"
)

# ------------------------------------------------------------------
# 1. Load isolate data (only the columns we need, to avoid the
#    crash caused by the column literally named "NA")
# ------------------------------------------------------------------
data <- read_excel("NGIsolates.xlsx", sheet = "NGIsolates") %>%
  select(IsolateID, FKSampleID, bact_genus, bact_species,
         ERYdiff, TYLdiff, VANdiff, CIPdiff, AMPdiff, TETdiff)

# Farm/timepoint metadata (used for the "before filtering" breakdown,
# since some Enterococcus samples might not be in the final valid-sample list)
sample_meta <- read_excel("NGSamples.xlsx") %>%
  select(FKSampleID, CalfID, FarmID, SamplingTimepoint)

# ------------------------------------------------------------------
# 2. Build the master valid-sample list from calves_per_farm_final.xlsx
#    (these are the final, faecal-only samples used in the analysis)
# ------------------------------------------------------------------
farm_sheets <- excel_sheets("calves_per_farm_final.xlsx")

valid_samples <- lapply(farm_sheets, function(sh) {
  df <- read_excel("calves_per_farm_final.xlsx", sheet = sh)
  df$CalfID <- as.character(df$CalfID)
  
  df %>%
    pivot_longer(cols = -CalfID, names_to = "SamplingLabel", values_to = "FKSampleID") %>%
    filter(!is.na(FKSampleID)) %>%
    mutate(
      FarmID_final = sh,
      SamplingTimepoint_final = as.integer(str_extract(SamplingLabel, "\\d+"))
    ) %>%
    select(FKSampleID, CalfNumber = CalfID, FarmID_final, SamplingTimepoint_final)
}) %>%
  bind_rows()

cat("Total valid (final, faecal) samples across all farms:", nrow(valid_samples), "\n\n")

# ------------------------------------------------------------------
# 3. Build the full Enterococcus faecalis/faecium isolate set
#    (BEFORE filtering), with Plate derived from IsolateID suffix
#    and Farm joined in
# ------------------------------------------------------------------
entero_all <- data %>%
  mutate(
    suffix = str_sub(IsolateID, -2, -1),
    Plate = case_when(
      suffix %in% c("ba", "bb") ~ "No selection",
      suffix %in% c("bc", "bd") ~ "Tylosin",
      suffix %in% c("be", "bf") ~ "Tetracycline",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(bact_genus == "Enterococcus",
         bact_species %in% c("faecalis", "faecium"),
         !is.na(Plate)) %>%
  left_join(sample_meta, by = "FKSampleID")

# ------------------------------------------------------------------
# 4. BEFORE-filtering counts, by farm and by plate
# ------------------------------------------------------------------
cat("===== BEFORE filtering =====\n")
cat("Total Enterococcus (faecalis/faecium) isolates:", nrow(entero_all), "\n")
cat("Total distinct samples:", n_distinct(entero_all$FKSampleID), "\n\n")

before_by_farm <- entero_all %>%
  distinct(FKSampleID, FarmID) %>%
  count(FarmID, name = "n_samples_before")

before_by_farm_isolates <- entero_all %>%
  count(FarmID, name = "n_isolates_before")

before_by_plate <- entero_all %>%
  distinct(FKSampleID, Plate) %>%
  count(Plate, name = "n_samples_before")

before_by_plate_isolates <- entero_all %>%
  count(Plate, name = "n_isolates_before")

print(before_by_farm)
print(before_by_farm_isolates)
print(before_by_plate)
print(before_by_plate_isolates)

# ------------------------------------------------------------------
# 5. Filter to ONLY the valid (final, faecal) samples
# ------------------------------------------------------------------
entero <- entero_all %>%
  semi_join(valid_samples, by = "FKSampleID")

entero$Plate <- factor(entero$Plate,
                       levels = c("No selection", "Tylosin", "Tetracycline"))

# ------------------------------------------------------------------
# 6. AFTER-filtering counts, and what was LOST, by farm and by plate
# ------------------------------------------------------------------
cat("\n===== AFTER filtering (valid/faecal samples only) =====\n")
cat("Total Enterococcus (faecalis/faecium) isolates:", nrow(entero), "\n")
cat("Total distinct samples:", n_distinct(entero$FKSampleID), "\n\n")

after_by_farm <- entero %>%
  distinct(FKSampleID, FarmID) %>%
  count(FarmID, name = "n_samples_after")

after_by_farm_isolates <- entero %>%
  count(FarmID, name = "n_isolates_after")

after_by_plate <- entero %>%
  distinct(FKSampleID, Plate) %>%
  count(Plate, name = "n_samples_after")

after_by_plate_isolates <- entero %>%
  count(Plate, name = "n_isolates_after")

lost_by_farm <- before_by_farm %>%
  full_join(after_by_farm, by = "FarmID") %>%
  mutate(across(starts_with("n_samples"), ~ replace_na(.x, 0)),
         n_samples_lost = n_samples_before - n_samples_after)

lost_by_farm_isolates <- before_by_farm_isolates %>%
  full_join(after_by_farm_isolates, by = "FarmID") %>%
  mutate(across(starts_with("n_isolates"), ~ replace_na(.x, 0)),
         n_isolates_lost = n_isolates_before - n_isolates_after)

lost_by_plate <- before_by_plate %>%
  full_join(after_by_plate, by = "Plate") %>%
  mutate(across(starts_with("n_samples"), ~ replace_na(.x, 0)),
         n_samples_lost = n_samples_before - n_samples_after)

lost_by_plate_isolates <- before_by_plate_isolates %>%
  full_join(after_by_plate_isolates, by = "Plate") %>%
  mutate(across(starts_with("n_isolates"), ~ replace_na(.x, 0)),
         n_isolates_lost = n_isolates_before - n_isolates_after)

cat("Samples lost by farm:\n"); print(lost_by_farm)
cat("\nIsolates lost by farm:\n"); print(lost_by_farm_isolates)
cat("\nSamples lost by plate:\n"); print(lost_by_plate)
cat("\nIsolates lost by plate:\n"); print(lost_by_plate_isolates)

write.csv(lost_by_farm, "enterococcus_samples_lost_by_farm.csv", row.names = FALSE)
write.csv(lost_by_farm_isolates, "enterococcus_isolates_lost_by_farm.csv", row.names = FALSE)
write.csv(lost_by_plate, "enterococcus_samples_lost_by_plate.csv", row.names = FALSE)
write.csv(lost_by_plate_isolates, "enterococcus_isolates_lost_by_plate.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 7. Resistance analysis - ONLY on the filtered (valid, faecal) set
#    ECOFF: faecium TYL <= 11, faecalis TYL <= 17 (species-specific).
#    All other antibiotics use the same cutoff for both species.
# ------------------------------------------------------------------
amr <- entero %>%
  mutate(
    ERY = ERYdiff, TYL = TYLdiff, VAN = VANdiff,
    CIP = CIPdiff, AMP = AMPdiff, TET = TETdiff
  ) %>%
  mutate(
    ERY_R = ifelse(ERY <= 13.01, 1, 0),
    TYL_R = case_when(
      bact_species == "faecium"  & TYL <= 11.01 ~ 1,
      bact_species == "faecium"  & TYL >  11.01 ~ 0,
      bact_species == "faecalis" & TYL <= 17.01 ~ 1,
      bact_species == "faecalis" & TYL >  17.01 ~ 0,
      TRUE ~ NA_real_
    ),
    VAN_R = ifelse(VAN <= 14.01, 1, 0),
    CIP_R = ifelse(CIP <= 15.01, 1, 0),
    AMP_R = ifelse(AMP <= 16.01, 1, 0),
    TET_R = ifelse(TET <= 14.01, 1, 0)
  )

resistance_cols <- names(amr)[str_ends(names(amr), "_R")]
antimicrobial_cols <- str_remove(resistance_cols, "_R$")

# ------------------------------------------------------------------
# 8. SAMPLE-level: which samples are positive for each antimicrobial,
#    split by farm and by species
# ------------------------------------------------------------------
sample_level <- amr %>%
  group_by(FKSampleID, FarmID, bact_species) %>%
  summarise(
    across(all_of(resistance_cols), ~ as.integer(any(.x == 1, na.rm = TRUE))),
    .groups = "drop"
  )

resistance_summary_samples <- sample_level %>%
  group_by(FarmID, bact_species) %>%
  summarise(
    n_samples = n(),
    across(all_of(resistance_cols), ~ sum(.x, na.rm = TRUE), .names = "{.col}_positive_samples"),
    .groups = "drop"
  )

print(resistance_summary_samples)
write.csv(resistance_summary_samples, "enterococcus_faecal_resistance_summary_samples.csv", row.names = FALSE)


# ------------------------------------------------------------------
# 9. ISOLATE-level UpSet plot of resistance profiles, split/coloured by farm
# ------------------------------------------------------------------
upset_data_isolate <- amr %>%
  select(IsolateID, FKSampleID, FarmID, bact_species, all_of(resistance_cols)) %>%
  rename_with(~ str_remove(.x, "_R$"), all_of(resistance_cols)) %>%
  mutate(across(all_of(antimicrobial_cols), ~ .x == 1))

# ---- E. faecalis ----
upset_faecalis <- upset_data_isolate %>% filter(bact_species == "faecalis")

p_upset_faecalis <- upset(
  upset_faecalis,
  antimicrobial_cols,
  name = "Resistant to",
  min_size = 1,
  width_ratio = 0.15,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = FarmID)
    ) + scale_fill_manual(values = farm_colors)
  )
) + plot_annotation(title = "E. faecalis resistance profile co-occurrence (isolate level, faecal samples only), coloured by farm")

print(p_upset_faecalis)
ggsave("faecalis_faecal_resistance_upset_isolate_by_farm.png", p_upset_faecalis,
       width = 12, height = 7, dpi = 300)

# ---- E. faecium ----
upset_faecium <- upset_data_isolate %>% filter(bact_species == "faecium")

p_upset_faecium <- upset(
  upset_faecium,
  antimicrobial_cols,
  name = "Resistant to",
  min_size = 1,
  width_ratio = 0.15,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = FarmID)
    ) + scale_fill_manual(values = farm_colors)
  )
) + plot_annotation(title = "E. faecium resistance profile co-occurrence (isolate level, faecal samples only), coloured by farm")

print(p_upset_faecium)
ggsave("faecium_faecal_resistance_upset_isolate_by_farm.png", p_upset_faecium,
       width = 12, height = 7, dpi = 300)

library(dplyr)

antimicrobial_cols <- c("ERY", "TYL", "VAN", "CIP", "AMP", "TET")

profiles <- upset_faecalis %>%
  rowwise() %>%
  mutate(
    Profile = paste(
      antimicrobial_cols[c_across(all_of(antimicrobial_cols))],
      collapse = "-"
    )
  ) %>%
  ungroup()

profiles %>%
  select(IsolateID, FarmID, Profile)
write_xlsx(profiles, "profiles.xlsx")
