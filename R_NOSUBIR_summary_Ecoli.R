#####################################################################
# Escherichia coli - AMR resistance analysis (isolate + sample level, milk handling)
#####################################################################

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(knitr)

setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/")

data <- read_excel("NGIsolates.xlsx")


# ------------------------------------------------------------------
# 1. Read every sheet (farm) and reshape wide -> long
# ------------------------------------------------------------------
farm_sheets <- excel_sheets("calves_per_farm_final.xlsx")

valid_samples <- lapply(farm_sheets, function(sh) {
  df <- read_excel("calves_per_farm_final.xlsx", sheet = sh)
  df$CalfID <- as.character(df$CalfID)   # some sheets have numeric CalfIDs, others text (e.g. "Black feet") - force consistent type before combining
  
  
  df %>%
    pivot_longer(
      cols = -CalfID,
      names_to = "SamplingLabel",
      values_to = "FKSampleID"
    ) %>%
    filter(!is.na(FKSampleID)) %>%
    mutate(
      FarmID = sh,
      SamplingTimepoint = as.integer(str_extract(SamplingLabel, "\\d+"))
    ) %>%
    select(FKSampleID, CalfNumber = CalfID, FarmID, SamplingTimepoint)
  # renamed to CalfNumber to avoid clashing with the CalfID column
  # in your other metadata file (that one encodes milk/pellet, not calf ID)
}) %>%
  bind_rows()

cat("Total valid samples across all farms:", nrow(valid_samples), "\n")
print(count(valid_samples, FarmID))

# Sanity check: any FKSampleID listed more than once? (data entry issue)
dupes <- valid_samples %>% filter(duplicated(FKSampleID) | duplicated(FKSampleID, fromLast = TRUE))
if (nrow(dupes) > 0) {
  cat("\nWARNING: duplicate FKSampleID entries found - check your spreadsheet:\n")
  print(dupes)
}

# ------------------------------------------------------------------
# 2. Restrict `data` to ONLY the isolates whose FKSampleID is in the
#    final valid-sample list. Do this right after loading data, before
#    any genus/species filtering, so every downstream plot/table
#    automatically respects it.
# ------------------------------------------------------------------
data <- data %>%
  semi_join(valid_samples, by = "FKSampleID")

cat("\nIsolates remaining after restricting to valid samples:", nrow(data), "\n")

sample_meta <- read_excel("NGSamples.xlsx") %>%
  select(FKSampleID, CalfID, FarmID, SamplingTimepoint)   # keep only what we need




# ------------------------------------------------------------------
# 1. Build amr: E. coli only, diffusion zones mapped in, metadata joined ONCE
# ------------------------------------------------------------------
amr <- data %>%
  filter(bact_genus == "Escherichia", bact_species == "coli") %>%
  mutate(
    AUG = AUGdiff,
    TS  = TSdiff,
    STR = STRdiff,
    CIP = CIPdiff,
    CAZ = CAZdiff,
    TET = TETdiff,
    CHL = CHLdiff,
    CTX = CTXdiff
  ) %>%
  left_join(sample_meta, by = "FKSampleID")

# ------------------------------------------------------------------
# 2. Resistance flags (no species split needed - only "coli" here)
# ------------------------------------------------------------------
amr <- amr %>%
  mutate(
    AUG_R = ifelse(AUG <= 13.01, 1, 0),
    CIP_R = ifelse(CIP <= 21.01, 1, 0),
    CAZ_R = ifelse(CAZ <= 17.01, 1, 0),
    CHL_R = ifelse(CHL <= 12.01, 1, 0),
    CTX_R = ifelse(CTX <= 22.01, 1, 0),
    STR_R = ifelse(STR <= 11.01, 1, 0),
    TS_R  = ifelse(TS  <= 10.01, 1, 0),
    TET_R = ifelse(TET <= 11.01, 1, 0)
  )

# ------------------------------------------------------------------
# 3. Pull out milk isolates FIRST (while amr still has them), with
#    Resistant/Not resistant flags - not raw diffusion zones
# ------------------------------------------------------------------
milk_isolates <- amr %>%
  filter(str_detect(str_to_lower(CalfID), "milk")) %>%
  mutate(
    AUG_result = ifelse(AUG_R == 1, "Resistant", "Not resistant"),
    CAZ_result = ifelse(CAZ_R == 1, "Resistant", "Not resistant"),
    CHL_result = ifelse(CHL_R == 1, "Resistant", "Not resistant"),
    CIP_result = ifelse(CIP_R == 1, "Resistant", "Not resistant"),
    CTX_result = ifelse(CTX_R == 1, "Resistant", "Not resistant"),
    STR_result = ifelse(STR_R == 1, "Resistant", "Not resistant"),
    TS_result  = ifelse(TS_R  == 1, "Resistant", "Not resistant"),
    TET_result = ifelse(TET_R == 1, "Resistant", "Not resistant")
  ) %>%
  select(IsolateID, FKSampleID, CalfID, FarmID, SamplingTimepoint,
         bact_genus, bact_species,
         AUG_result, CAZ_result, CHL_result, CIP_result, CTX_result,
         STR_result, TS_result, TET_result)

print(milk_isolates)
write.csv(milk_isolates, "ecoli_milk_isolates_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 4. NOW exclude milk from amr for the main (non-milk) analysis
# ------------------------------------------------------------------
amr <- amr %>%
  filter(!str_detect(str_to_lower(CalfID), "milk"))

# ------------------------------------------------------------------
# 5. Isolate-level resistance summary
# ------------------------------------------------------------------
resistance_summary <- amr %>%
  group_by(bact_species) %>%
  summarise(
    n_isolates    = n(),
    AUG_resistant = sum(AUG_R, na.rm = TRUE),
    CAZ_resistant = sum(CAZ_R, na.rm = TRUE),
    CHL_resistant = sum(CHL_R, na.rm = TRUE),
    CIP_resistant = sum(CIP_R, na.rm = TRUE),
    CTX_resistant = sum(CTX_R, na.rm = TRUE),
    STR_resistant = sum(STR_R, na.rm = TRUE),
    TS_resistant  = sum(TS_R, na.rm = TRUE),
    TET_resistant = sum(TET_R, na.rm = TRUE)
  )

print(resistance_summary)
write.csv(resistance_summary, "ecoli_resistance_summary.csv", row.names = FALSE)

# Long format, isolate level - all 8 antibiotics included
resistance_long <- amr %>%
  select(bact_species, AUG_R, CAZ_R, CHL_R, CIP_R, CTX_R, STR_R, TS_R, TET_R) %>%
  pivot_longer(cols = -bact_species, names_to = "Antimicrobial", values_to = "Resistant") %>%
  mutate(Antimicrobial = str_remove(Antimicrobial, "_R")) %>%
  group_by(bact_species, Antimicrobial) %>%
  summarise(
    n_tested      = sum(!is.na(Resistant)),
    n_resistant   = sum(Resistant, na.rm = TRUE),
    pct_resistant = 100 * n_resistant / n_tested,
    .groups = "drop"
  )

p_resistance_isolate <- ggplot(resistance_long, aes(x = Antimicrobial, y = pct_resistant)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(round(pct_resistant, 1), "%\n(", n_resistant, "/", n_tested, ")")),
            vjust = -0.1, size = 3) +
  labs(
    title = "Percentage of resistant E. coli isolates by antimicrobial",
    x = "Antimicrobial", y = "% resistant"
  ) +
  theme_bw()

print(p_resistance_isolate)
ggsave("ecoli_resistance_by_antimicrobial_isolate_level.png", p_resistance_isolate,
       width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 6. Sample-level resistance (a sample counts once per drug if ANY
#    of its isolates were resistant to that drug)
# ------------------------------------------------------------------
sample_level <- amr %>%
  group_by(FKSampleID, bact_species) %>%
  summarise(
    n_isolates = n(),
    AUG_R = as.integer(any(AUG_R == 1, na.rm = TRUE)),
    CAZ_R = as.integer(any(CAZ_R == 1, na.rm = TRUE)),
    CHL_R = as.integer(any(CHL_R == 1, na.rm = TRUE)),
    CIP_R = as.integer(any(CIP_R == 1, na.rm = TRUE)),
    CTX_R = as.integer(any(CTX_R == 1, na.rm = TRUE)),
    STR_R = as.integer(any(STR_R == 1, na.rm = TRUE)),
    TS_R  = as.integer(any(TS_R  == 1, na.rm = TRUE)),
    TET_R = as.integer(any(TET_R == 1, na.rm = TRUE)),
    .groups = "drop"
  )

resistance_summary_samples <- sample_level %>%
  group_by(bact_species) %>%
  summarise(
    n_samples             = n(),
    AUG_resistant_samples = sum(AUG_R, na.rm = TRUE),
    CAZ_resistant_samples = sum(CAZ_R, na.rm = TRUE),
    CHL_resistant_samples = sum(CHL_R, na.rm = TRUE),
    CIP_resistant_samples = sum(CIP_R, na.rm = TRUE),
    CTX_resistant_samples = sum(CTX_R, na.rm = TRUE),
    STR_resistant_samples = sum(STR_R, na.rm = TRUE),
    TS_resistant_samples  = sum(TS_R, na.rm = TRUE),
    TET_resistant_samples = sum(TET_R, na.rm = TRUE)
  )

print(resistance_summary_samples)
write.csv(resistance_summary_samples, "ecoli_resistance_summary_samples.csv", row.names = FALSE)

resistance_long_samples <- sample_level %>%
  select(bact_species, AUG_R, CAZ_R, CHL_R, CIP_R, CTX_R, STR_R, TS_R, TET_R) %>%
  pivot_longer(cols = -bact_species, names_to = "Antimicrobial", values_to = "Resistant") %>%
  mutate(Antimicrobial = str_remove(Antimicrobial, "_R")) %>%
  group_by(bact_species, Antimicrobial) %>%
  summarise(
    n_samples_tested    = sum(!is.na(Resistant)),
    n_resistant_samples = sum(Resistant, na.rm = TRUE),
    pct_resistant       = 100 * n_resistant_samples / n_samples_tested,
    .groups = "drop"
  )

p_resistance_sample <- ggplot(resistance_long_samples, aes(x = Antimicrobial, y = n_resistant_samples)) +
  geom_col(fill = "darkorange") +
  geom_text(aes(label = paste0(n_resistant_samples, "/", n_samples_tested)),
            vjust = -0.3, size = 3.5) +
  labs(
    title = "Number of E. coli samples resistant by antimicrobial",
    x = "Antimicrobial", y = "Number of resistant samples"
  ) +
  theme_bw()

print(p_resistance_sample)
ggsave("ecoli_resistance_by_antimicrobial_sample_level.png", p_resistance_sample,
       width = 8, height = 5, dpi = 300)
