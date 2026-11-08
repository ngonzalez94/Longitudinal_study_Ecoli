# ============================================================
# Antimicrobial Resistance Analysis - Enterococcus
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ---- 1. Load Data ---- 
data <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/Ecoli_Complete.xlsx", sheet=1)

# ---- 2. Assign Plate Group from IsolateID suffix ----
data <- data %>%
  mutate(
    plate_code = substr(IsolateID, nchar(IsolateID) - 1, nchar(IsolateID)),
    PlateGroup = case_when(
      plate_code %in% c("aa", "ab", "FE", "b'", "a'") ~ "No Antibiotic",
      plate_code %in% c("ac", "ad") ~ "Ceftazidime",
      plate_code %in% c("ae", "af", "fe") ~ "Cefotaxime",
      plate_code %in% c("ag", "ah") ~ "Tetracycline",
      TRUE ~ "Unknown"
    ),
    PlateGroup = factor(PlateGroup, levels = c("No Antibiotic", "Ceftazidime", "Cefotaxime", "Tetracycline"))
  )

# ---- 3. Apply Resistance Cut-offs (1 = Resistant, 0 = Susceptible) ----
ab_cols  <- c("AUG_R", "CIP_R", "CAZ_R", "TET_R", "CTX_R", "CHL_R", "STR_R", "TS_R")
ab_labels <- c("Ampicillin", "Ciprofloxacin", "Ceftazidime","Tetracycline", "Cefotaxime", "Chloramphenicol", "Streptomycin", "Trim/Sulf")


data <- data %>%
  mutate(
    AUG_R = ifelse(AUG <= 13.01, 1, 0),
    CAZ_R = ifelse(CAZ <= 17.01, 1, 0),
    CHL_R = ifelse(CHL <= 12.01, 1, 0),
    CIP_R = ifelse(CIP <= 21.01, 1, 0),
    CTX_R = ifelse(CTX <= 22.01, 1, 0),
    STR_R = ifelse(STR <= 11.01, 1, 0),
    TET_R = ifelse(TET <= 11.01, 1, 0),
    TS_R = ifelse(TS <= 10.01, 1, 0),
    # Resistant to at least one antibiotic (TYL11 faecium, TYL17 faecalis)
    any_R = ifelse(rowSums(across(all_of(ab_cols)), na.rm = TRUE) >= 1, 1, 0)
  )



# ============================================================
# REPORT SECTION 1: Isolate counts per plate
# ============================================================
cat("\n========================================\n")
cat("ISOLATE COUNTS PER PLATE\n")
cat("========================================\n")

plate_counts <- data %>%
  count(PlateGroup, name = "n_isolates") %>%
  mutate(Pct = round(n_isolates / sum(n_isolates) * 100, 1))

print(plate_counts)

# ============================================================
# REPORT SECTION 2: Resistant to at least one antibiotic per plate
# ============================================================
cat("\n========================================\n")
cat("RESISTANCE TO AT LEAST ONE ANTIBIOTIC\n")
cat("========================================\n")

any_resist_summary <- data %>%
  group_by(PlateGroup) %>%
  summarise(
    n_isolates   = n(),
    n_resistant  = sum(any_R, na.rm = TRUE),
    pct_resistant = round(n_resistant / n_isolates * 100, 1),
    .groups = "drop"
  )

print(any_resist_summary)

cat("\nFisher's exact test (any resistance ~ plate group):\n")
print(fisher.test(table(data$PlateGroup, data$any_R), simulate.p.value = TRUE))

# Use Fisher's exact if any expected count < 5
if (any(chisq_any$expected < 5)) {
  cat("\nNote: Expected counts < 5 detected. Using Fisher's exact test instead:\n")
  print(fisher.test(table(data$PlateGroup, data$any_R), simulate.p.value = TRUE))
}

# ============================================================
# REPORT SECTION 3: Per-antibiotic resistance by plate group
# ============================================================
cat("\n========================================\n")
cat("RESISTANCE BY ANTIBIOTIC AND PLATE GROUP\n")
cat("========================================\n")

per_ab_summary <- data %>%
  group_by(PlateGroup) %>%
  summarise(
    n_isolates = n(),
    across(all_of(ab_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = all_of(ab_cols),
               names_to = "Antibiotic", values_to = "n_Resistant") %>%
  mutate(
    Antibiotic = factor(Antibiotic, levels = ab_cols, labels = ab_labels),
    Pct = round(n_Resistant / n_isolates * 100, 1)
  )

print(per_ab_summary %>% select(PlateGroup, Antibiotic, n_isolates, n_Resistant, Pct))

# ============================================================
# REPORT SECTION 4: Fisher's exact test per antibiotic
# ============================================================
cat("\n========================================\n")
cat("STATISTICAL TESTS PER ANTIBIOTIC (Fisher's Exact)\n")
cat("========================================\n")

fisher_results <- lapply(seq_along(ab_cols), function(i) {
  ab  <- ab_cols[i]
  lab <- ab_labels[i]
  tbl <- table(data$PlateGroup, data[[ab]])
  
  # Skip if table doesn't have at least 2 rows and 2 columns
  if (nrow(tbl) < 2 || ncol(tbl) < 2) {
    return(data.frame(
      Antibiotic  = lab,
      p_value     = NA,
      Significant = "not testable (no variation)"
    ))
  }
  
  res <- fisher.test(tbl, simulate.p.value = TRUE)
  
  data.frame(
    Antibiotic  = lab,
    p_value     = round(res$p.value, 4),
    Significant = ifelse(res$p.value < 0.05, "YES *", "no")
  )
})

chisq_df <- bind_rows(fisher_results) %>%
  arrange(p_value)

print(chisq_df)

#find antibiotic only one column 
for (ab in ab_cols) {
  cat(ab, "-> unique values:", unique(data[[ab]]), "\n")
}
# ============================================================
# PLOTS
# ============================================================
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/")
# Plot 1: Isolate counts per plate
p1 <- ggplot(plate_counts, aes(x = PlateGroup, y = n_isolates, fill = PlateGroup)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(n_isolates, "\n(", Pct, "%)")), 
            vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("No Antibiotic" = "darkgreen",
                               "Ceftazidime"       = "lightblue",
                               "Cefotaxime"       = "darkblue",
                               "Tetracycline"  = "darkred")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Number of Isolates per Selection Plate",
       x = "Plate Group", y = "Number of Isolates") +
  theme_bw(base_size = 13)

print(p1)
ggsave("plot1_isolate_counts_ecoli.png", p1, width = 7, height = 5, dpi = 300)

# Plot 2: % resistant to at least one antibiotic per plate
p2 <- ggplot(any_resist_summary, 
             aes(x = PlateGroup, y = pct_resistant, fill = PlateGroup)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(pct_resistant, "%\n(", n_resistant, "/", n_isolates, ")")),
            vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("No Antibiotic" = "darkgreen",
                               "Ceftazidime"       = "lightblue",
                               "Cefotaxime"       = "darkblue",
                               "Tetracycline"  = "darkred")) +
  scale_y_continuous(limits = c(0, 110), labels = label_percent(scale = 1)) +
  labs(title = "Isolates Resistant to ≥1 Antibiotic by Selection Plate",
       x = "Plate Group", y = "% Resistant") +
  theme_bw(base_size = 13)

print(p2)
ggsave("plot2_any_resistance_ecoli.png", p2, width = 7, height = 5, dpi = 300)

# Plot 3: Per-antibiotic resistance % by plate group (grouped bar)
p3 <- ggplot(per_ab_summary, 
             aes(x = Antibiotic, y = Pct, fill = PlateGroup)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0(Pct, "%")),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("No Antibiotic" = "darkgreen",
                               "Ceftazidime"       = "lightblue",
                               "Cefotaxime"       = "darkblue",
                               "Tetracycline"  = "darkred")) +
  scale_y_continuous(limits = c(0, 110), labels = label_percent(scale = 1)) +
  labs(#title    = "Resistance per Antibiotic by Selection Plate",
       #subtitle = "* = significant difference between plates (p < 0.05)",
       x = "Antibiotic", y = "% Resistant", fill = "Plate Group") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Add significance stars from chisq_df
sig_abs <- chisq_df %>% filter(Significant == "YES *") %>% pull(Antibiotic)
if (length(sig_abs) > 0) {
  sig_data <- per_ab_summary %>%
    filter(Antibiotic %in% sig_abs) %>%
    group_by(Antibiotic) %>%
    summarise(y_pos = max(Pct) + 8, .groups = "drop")
  
  p3 <- p3 + geom_text(data = sig_data, 
                       aes(x = Antibiotic, y = y_pos, label = "*"),
                       inherit.aes = FALSE, size = 6, color = "red")
}

print(p3)
ggsave("plot3_per_antibiotic_by_plate_ecoli.png", p3, width = 9, height = 5, dpi = 300)
