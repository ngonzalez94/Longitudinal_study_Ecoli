# ==============================================================================
# Virulence Gene Heatmap
# Sources: Abricate (yellow), VirulenceFinder (blue), Both (green)
# Layout:  per-farm with Diet Stage row splits + combined all-farms
# Metadata from ARG_and_Metadata.xlsx (Farm, Diet_Stage, ST, Resistance_Profile)
# ==============================================================================

library(ComplexHeatmap)
library(circlize)
library(readxl)
library(dplyr)
library(tidyr)
library(grid)

# ==============================================================================
# 1. FILE PATHS — UPDATE THESE
# ==============================================================================

metadata_path      <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/ARG_and_Metadata.xlsx"
metadata_sheet     <- 3          # sheet containing isolate metadata

# Abricate virulence summary file (wide format: isolates as rows, genes as columns)
# Same structure as ARG_and_Metadata.xlsx — #FILE column + gene columns with % or "."
abricate_vir_path  <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/summary_vfdb.xlsx"

# VirulenceFinder file (long format: one row per gene hit)
# Expected columns: Sample, Virulence factor
virfinder_path     <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/virulence_finder.xlsx"

save_files <- TRUE

# ==============================================================================
# 2. COLORS
# ==============================================================================

# Cell colors: 0=absent, 1=Abricate only, 2=VirulenceFinder only, 3=Both
source_colors <- colorRamp2(c(0, 1, 2, 3),
                            c("white", "#FFD700", "#2166AC", "#33A02C"))

# Resistance Profile label mapping
rp_label_map <- c(
  "A" = "TET-TS-STR",
  "B" = "TET-TS-STR-CTX",
  "C" = "TET-TS-AUG",
  "D" = "TET-TS-AUG-CTX",
  "E" = "TET-STR-CHL",
  "F" = "TET-CTX-STR",
  "G" = "TET-AUG-STR"
)

custom_rp_colors <- c(
  "TET-TS-STR"     = "#e41a1c",
  "TET-TS-STR-CTX" = "#377eb8",
  "TET-TS-AUG"     = "#4daf4a",
  "TET-TS-AUG-CTX" = "#984ea3",
  "TET-STR-CHL"    = "#ff7f00",
  "TET-CTX-STR"    = "#ffff33",
  "TET-AUG-STR"    = "#a65628"
)

custom_st_colors <- c(
  "10"   = "#a50026", "1613" = "#f46d43", "227"  = "#fee08b",
  "23"   = "#ffffbf", "3572" = "#d9ef8b", "362"  = "#66bd63",
  "48"   = "#006837", "58"   = "#053061", "69"   = "#2166ac",
  "783"  = "#4393c3", "8102" = "#92c5de", "85"   = "#d1e5f0",
  "88"   = "#fde0ef", "4559" = "#de77ae", "56"   = "#c51b7d",
  "117"  = "#8e0152", "NM"   = "grey"
)

diet_colors <- c(
  "Pre-WM"  = "#ADD8E6",
  "WM"      = "#1F78B4",
  "Post-WM" = "#B2E0A0"
)

# ==============================================================================
# 3. LOAD & PREPARE METADATA
# ==============================================================================

meta_raw <- read_excel(metadata_path, sheet = metadata_sheet)

# Auto-detect isolate column
candidate_cols <- c("#FILE", "Isolate", "Sample", "isolate", "sample",
                    "Strain", "strain", "ID")
isolate_col <- NULL
for (cand in candidate_cols) {
  if (cand %in% colnames(meta_raw)) { isolate_col <- cand; break }
}
if (is.null(isolate_col)) {
  cat("ERROR: Cannot find isolate column. Columns available:\n")
  print(colnames(meta_raw))
  stop("Set isolate_col manually.")
}
cat("Isolate column detected: '", isolate_col, "'\n", sep = "")

isolates <- gsub("_contig\\.fa$", "", as.character(meta_raw[[isolate_col]]))

metadata <- as.data.frame(meta_raw[, c("Resistance_Profile", "Farm",
                                       "Diet_Stage", "ST")])
rownames(metadata) <- isolates

# Translate resistance profile codes to labels
metadata$Resistance_Profile <- rp_label_map[as.character(metadata$Resistance_Profile)]
metadata$Resistance_Profile[is.na(metadata$Resistance_Profile)] <-
  as.character(meta_raw$Resistance_Profile[is.na(metadata$Resistance_Profile)])

# Set Diet Stage order
metadata$Diet_Stage <- factor(as.character(metadata$Diet_Stage),
                              levels = c("Pre-WM", "WM", "Post-WM"),
                              ordered = TRUE)

# Metadata color schemes
unique_profiles <- sort(unique(metadata$Resistance_Profile))
profile_colors  <- custom_rp_colors
missing_rp <- setdiff(unique_profiles, names(custom_rp_colors))
if (length(missing_rp) > 0)
  profile_colors <- c(profile_colors,
                      setNames(hcl.colors(length(missing_rp), "Set 2"), missing_rp))

unique_farms <- sort(unique(metadata$Farm))
farm_colors  <- setNames(rainbow(length(unique_farms)), unique_farms)

unique_st <- sort(unique(as.character(metadata$ST)))
st_colors  <- custom_st_colors
missing_st <- setdiff(unique_st, names(custom_st_colors))
if (length(missing_st) > 0)
  st_colors <- c(st_colors,
                 setNames(hcl.colors(length(missing_st), "Pastel 1"), missing_st))

all_isolates_df <- data.frame(Sample = isolates, stringsAsFactors = FALSE)

# ==============================================================================
# 4. LOAD ABRICATE VIRULENCE DATA (wide/summary format)
#    Structure: #FILE column + one column per gene, values = % identity or "."
# ==============================================================================

ab_vir_raw <- read_excel(abricate_vir_path)

# Clean isolate names from #FILE column
ab_vir_raw[["#FILE"]] <- gsub("_contig\\.fa$", "",
                              as.character(ab_vir_raw[["#FILE"]]))

# All columns except #FILE and NUM_FOUND are gene columns
ab_gene_cols  <- setdiff(colnames(ab_vir_raw), c("#FILE", "NUM_FOUND"))
ab_gene_data  <- as.data.frame(ab_vir_raw[, ab_gene_cols])

# Presence/absence: "." or NA = 0, any value = 1
ab_vir_pa <- apply(ab_gene_data, 2, function(x)
  ifelse(x == "." | is.na(x), 0L, 1L))

ab_vir_wide   <- data.frame(Sample = ab_vir_raw[["#FILE"]],
                            ab_vir_pa, stringsAsFactors = FALSE)
ab_vir_wide   <- all_isolates_df %>% left_join(ab_vir_wide, by = "Sample")
ab_vir_matrix <- as.matrix(ab_vir_wide[, -1, drop = FALSE])
ab_vir_matrix[is.na(ab_vir_matrix)] <- 0L
rownames(ab_vir_matrix) <- isolates

# ==============================================================================
# 5. LOAD VIRULENCEFINDER DATA
# ==============================================================================

vf_raw <- read_excel(virfinder_path)

# VirulenceFinder uses Sample and "Virulence factor" columns
vf_raw$Sample <- gsub("_contig\\.fa$", "", as.character(vf_raw$Sample))

vf_long <- vf_raw %>%
  select(Sample, Gene = `Virulence factor`) %>%
  distinct() %>%
  mutate(present = 1L)

vf_wide   <- vf_long %>%
  pivot_wider(names_from = Gene, values_from = present, values_fill = 0L)

vf_wide   <- all_isolates_df %>% left_join(vf_wide, by = "Sample")
vf_matrix <- as.matrix(vf_wide[, -1, drop = FALSE])
vf_matrix[is.na(vf_matrix)] <- 0L
rownames(vf_matrix) <- isolates

cat("=== Virulence Diagnostics ===\n")
cat("Abricate virulence genes  :", ncol(ab_vir_matrix), "\n")
cat("VirulenceFinder genes     :", ncol(vf_matrix), "\n")

# ==============================================================================
# 6. BUILD COMBINED SOURCE MATRIX
#    0=absent | 1=Abricate only | 2=VirulenceFinder only | 3=Both
# ==============================================================================

all_vir_genes <- union(colnames(ab_vir_matrix), colnames(vf_matrix))

expand_mat <- function(mat, all_genes, row_names) {
  out <- matrix(0L, nrow = length(row_names), ncol = length(all_genes),
                dimnames = list(row_names, all_genes))
  shared <- intersect(colnames(mat), all_genes)
  out[, shared] <- mat[, shared]
  out
}

ab_full <- expand_mat(ab_vir_matrix, all_vir_genes, isolates)
vf_full <- expand_mat(vf_matrix,     all_vir_genes, isolates)

vir_source_matrix <- ab_full + 2L * vf_full

# Sort genes alphabetically (no category splitting)
vir_source_matrix <- vir_source_matrix[, sort(colnames(vir_source_matrix)),
                                       drop = FALSE]

cat("Total unique virulence genes:", ncol(vir_source_matrix), "\n\n")

# ==============================================================================
# 7. SOURCE LEGEND
# ==============================================================================

source_legend <- Legend(
  title     = "Virulence Detection Source",
  labels    = c("Absent", "Abricate", "VirulenceFinder", "Both"),
  legend_gp = gpar(fill = c("white", "#FFD700", "#2166AC", "#33A02C"),
                   col  = "gray60"),
  border    = "gray60"
)

# ==============================================================================
# 8. HELPER: ROW ANNOTATION
# ==============================================================================

make_row_ha <- function(meta, include_farm = FALSE) {
  if (include_farm) {
    rowAnnotation(
      Farm               = meta$Farm,
      `Resistance Profile` = meta$Resistance_Profile,
      `Diet Stage`       = meta$Diet_Stage,
      ST                 = as.character(meta$ST),
      col = list(Farm               = farm_colors,
                 `Resistance Profile` = profile_colors,
                 `Diet Stage`       = diet_colors,
                 ST                 = st_colors),
      annotation_name_side = "top",
      annotation_name_gp   = gpar(fontsize = 10),
      simple_anno_size     = unit(0.5, "cm"),
      gap                  = unit(1, "mm")
    )
  } else {
    rowAnnotation(
      `Resistance Profile` = meta$Resistance_Profile,
      `Diet Stage`       = meta$Diet_Stage,
      ST                 = as.character(meta$ST),
      col = list(`Resistance Profile` = profile_colors,
                 `Diet Stage`       = diet_colors,
                 ST                 = st_colors),
      annotation_name_side = "top",
      annotation_name_gp   = gpar(fontsize = 10),
      simple_anno_size     = unit(0.5, "cm"),
      gap                  = unit(1, "mm")
    )
  }
}

# ==============================================================================
# 9. PER-FARM HEATMAPS
# ==============================================================================

farms <- sort(unique(metadata$Farm))
cat("=== Creating per-farm virulence heatmaps ===\n")

for (farm in farms) {
  
  cat("\nProcessing Farm", farm, "...\n")
  
  idx        <- which(metadata$Farm == farm)
  farm_meta  <- metadata[idx, , drop = FALSE]
  farm_mat   <- vir_source_matrix[idx, , drop = FALSE]
  
  if (nrow(farm_mat) == 0) { cat("  No isolates — skipping.\n"); next }
  
  # Filter to genes present in this farm
  keep      <- colSums(farm_mat > 0) > 0
  farm_mat  <- farm_mat[, keep, drop = FALSE]
  cat("  Virulence genes present:", sum(keep), "\n")
  
  # Re-apply ordered Diet Stage factor
  farm_meta$Diet_Stage <- factor(as.character(farm_meta$Diet_Stage),
                                 levels = c("Pre-WM", "WM", "Post-WM"),
                                 ordered = TRUE)
  
  row_ha <- make_row_ha(farm_meta, include_farm = FALSE)
  
  ht <- Heatmap(
    farm_mat,
    name                = "Detection Source",
    col                 = source_colors,
    show_heatmap_legend = FALSE,
    
    row_names_side  = "left",
    row_names_gp    = gpar(fontsize = 8),
    row_title       = paste("Farm", farm, "- Isolates"),
    row_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows    = TRUE,
    show_row_dend   = TRUE,
    row_dend_width  = unit(2, "cm"),
    row_split       = factor(as.character(farm_meta$Diet_Stage),
                             levels = c("Pre-WM", "WM", "Post-WM"),
                             ordered = TRUE),
    row_gap         = unit(1, "mm"),
    
    column_names_side = "bottom",
    column_names_gp   = gpar(fontsize = 7),
    column_names_rot  = 45,
    column_title      = "Virulence Genes",
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    cluster_columns   = TRUE,
    show_column_dend  = TRUE,
    column_dend_height = unit(2, "cm"),
    
    left_annotation = row_ha,
    border          = TRUE,
    rect_gp         = gpar(col = "gray90", lwd = 0.5)
  )
  
  draw(ht,
       annotation_legend_side = "right",
       merge_legend           = TRUE,
       annotation_legend_list = list(source_legend))
  
  if (save_files) {
    pdf_out <- paste0("virulence_heatmap_Farm_", farm, ".pdf")
    png_out <- paste0("virulence_heatmap_Farm_", farm, ".png")
    pdf(pdf_out, width = 20, height = 10)
    draw(ht, annotation_legend_side = "right", merge_legend = TRUE,
         annotation_legend_list = list(source_legend))
    dev.off()
    png(png_out, width = 2000, height = 1000, res = 100)
    draw(ht, annotation_legend_side = "right", merge_legend = TRUE,
         annotation_legend_list = list(source_legend))
    dev.off()
    cat("  Saved:", pdf_out, "&", png_out, "\n")
  }
}

# ==============================================================================
# 10. COMBINED ALL-FARMS HEATMAP
# ==============================================================================

cat("\n=== Creating combined all-farms virulence heatmap ===\n")

keep_all  <- colSums(vir_source_matrix > 0) > 0
comb_mat  <- vir_source_matrix[, keep_all, drop = FALSE]
cat("Virulence genes across all farms:", sum(keep_all), "\n")

metadata$Farm <- factor(metadata$Farm, levels = sort(unique(metadata$Farm)))

row_ha_all <- make_row_ha(metadata, include_farm = TRUE)

ht_all <- Heatmap(
  comb_mat,
  name                = "Detection Source",
  col                 = source_colors,
  show_heatmap_legend = FALSE,
  
  row_names_side     = "left",
  row_names_gp       = gpar(fontsize = 7),
  row_title          = "Isolates (All Farms)",
  row_title_gp       = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows       = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend      = TRUE,
  row_dend_width     = unit(2, "cm"),
  row_split          = metadata$Farm,
  row_gap            = unit(3, "mm"),
  
  column_names_side  = "bottom",
  column_names_gp    = gpar(fontsize = 7),
  column_names_rot   = 45,
  column_title       = "Virulence Genes",
  column_title_gp    = gpar(fontsize = 12, fontface = "bold"),
  cluster_columns    = TRUE,
  show_column_dend   = TRUE,
  column_dend_height = unit(2, "cm"),
  
  left_annotation    = row_ha_all,
  border             = TRUE,
  rect_gp            = gpar(col = "gray90", lwd = 0.5)
)

draw(ht_all,
     annotation_legend_side = "right",
     merge_legend           = TRUE,
     annotation_legend_list = list(source_legend))

if (save_files) {
  pdf("virulence_heatmap_ALL_FARMS.pdf", width = 24, height = 16)
  draw(ht_all, annotation_legend_side = "right", merge_legend = TRUE,
       annotation_legend_list = list(source_legend))
  dev.off()
  png("virulence_heatmap_ALL_FARMS.png", width = 2400, height = 1600, res = 100)
  draw(ht_all, annotation_legend_side = "right", merge_legend = TRUE,
       annotation_legend_list = list(source_legend))
  dev.off()
  cat("Saved: virulence_heatmap_ALL_FARMS.pdf & virulence_heatmap_ALL_FARMS.png\n")
}

cat("\n=== ALL DONE! ===\n")