# ==============================================================================
# Antimicrobial Resistance Gene Heatmap with Metadata Annotations
# ==============================================================================
# This script creates a heatmap showing gene presence/absence across isolates,
# with genes organized by antimicrobial class and isolates annotated with
# farm, resistance profile, and diet stage metadata.
# ==============================================================================

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(dplyr)
library(tidyr)
library(grid)  # Needed for drawing circles

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

#sheet1 resfinder, sheet 2 card, sheet 3 ncbi
main_data <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/ARG_and_Metadata.xlsx", sheet = 3)

gene_classes <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/Ecoli_mechanisms.xlsx")

# =============================================================================
# CUSTOMIZE COLORS HERE - Edit these to change annotation colors
# =============================================================================

# REPLACE THE NULL VALUES BELOW WITH YOUR COLOR DEFINITIONS

# Resistance Profile colors - EDIT THIS:
#custom_rp_colors <- NULL  # ← REPLACE this NULL with your color vector

# Example (uncomment and modify):
custom_rp_colors <- c(
  "A" = "#e41a1c",
  "B" = "#377eb8", 
  "C" = "#4daf4a",
  "D" = "#984ea3", 
  "E" = "#ff7f00", 
  "F" = "#ffff33",
  "G" = "#a65628"
)

# Sequence Type colors - EDIT THIS:
#custom_st_colors <- NULL  # ← REPLACE this NULL with your color vector

# ST colors
custom_st_colors <- c(
  "10" = "#a50026",
  "1613" = "#f46d43",
  "227" = "#fee08b",
  "23" = "#ffffbf",
  "3572" = "#d9ef8b", 
  "362" = "#66bd63", 
  "48" = "#006837", 
  "58" = "#053061", 
  "69" = "#2166ac", 
  "783" = "#4393c3", 
  "8102" = "#92c5de", 
  "85" = "#d1e5f0", 
  "88" = "#fde0ef", 
  "4559" = "#de77ae",
  "56" = "#c51b7d", 
  "117" = "#8e0152",
  "NM" = "grey"
)

# =============================================================================

# SAVE FILES OR JUST PREVIEW?
save_files <- TRUE  # Set to FALSE to just preview plots without saving files

# =============================================================================
# ==============================================================================
# 2. PREPARE DATA
# ==============================================================================

# Extract metadata columns
metadata_cols <- c("Resistance_Profile", "Farm", "Diet_Stage", "ST")
isolate_col <- "Isolate"

# Extract isolate names
isolates <- main_data[[isolate_col]]

# Clean isolate names - remove "_contig.fa" suffix
isolates <- gsub("_contig\\.fa$", "", isolates)

# Extract metadata
metadata <- main_data[, metadata_cols]
rownames(metadata) <- isolates

# Extract gene data (all columns except isolate, NUM_FOUND, and metadata)
exclude_cols <- c(isolate_col, "NUM_FOUND", metadata_cols)
gene_cols <- setdiff(colnames(main_data), exclude_cols)
gene_data <- main_data[, gene_cols]

# Convert gene data to presence/absence matrix
# "." means absent, any numeric value means present
gene_matrix <- apply(gene_data, 2, function(x) {
  ifelse(x == "." | is.na(x), 0, 1)
})
rownames(gene_matrix) <- isolates

# Optional: If you want to use the actual percentage values instead of 0/1:
# Uncomment the lines below and comment out the above gene_matrix code
# gene_matrix <- apply(gene_data, 2, function(x) {
#   ifelse(x == "." | is.na(x), 0, as.numeric(x))
# })
# rownames(gene_matrix) <- isolates

# ==============================================================================
# 3. ORGANIZE GENES BY ANTIMICROBIAL CLASS
# ==============================================================================

# Match genes in your data with the gene classes
# Create a data frame with gene names and their classes
gene_info <- data.frame(
  Genes = gene_cols,
  stringsAsFactors = FALSE
)

# Merge with gene classes (using "Genes" and "Classification" columns)
gene_info <- left_join(gene_info, gene_classes, by = "Genes")

# Handle any genes without class assignment
gene_info$Classification[is.na(gene_info$Classification)] <- "Unknown"

# Order genes by antimicrobial class
gene_info <- gene_info[order(gene_info$Classification), ]

# Reorder gene matrix columns to match
gene_matrix <- gene_matrix[, gene_info$Genes]

# ==============================================================================
# 4. CREATE COLOR SCHEMES FOR METADATA
# ==============================================================================

# Define colors for Resistance Profile
unique_profiles <- sort(unique(metadata$Resistance_Profile))
if (!is.null(custom_rp_colors)) {
  # Use custom colors if provided, but fill in missing ones
  profile_colors <- custom_rp_colors
  # Add auto colors for any missing profiles
  missing_profiles <- setdiff(unique_profiles, names(custom_rp_colors))
  if (length(missing_profiles) > 0) {
    cat("WARNING: Adding auto colors for missing Resistance Profiles:", paste(missing_profiles, collapse=", "), "\n")
    auto_colors <- setNames(
      hcl.colors(length(missing_profiles), palette = "Set 2"),
      missing_profiles
    )
    profile_colors <- c(profile_colors, auto_colors)
  }
  cat("Using custom colors for Resistance Profile\n")
} else {
  # Auto-generate colors for all profiles
  profile_colors <- setNames(
    hcl.colors(length(unique_profiles), palette = "Set 2"),
    unique_profiles
  )
}

# Define colors for Farm
unique_farms <- unique(metadata$Farm)
farm_colors <- setNames(
  rainbow(length(unique_farms)),
  unique_farms
)

# Define colors for Diet Stage
unique_diets <- sort(unique(metadata$Diet_Stage))
diet_colors <- setNames(
  hcl.colors(length(unique_diets), palette = "Dark 3"),
  unique_diets
)

# Define colors for ST (Sequence Type)
unique_st <- sort(unique(metadata$ST))
if (!is.null(custom_st_colors)) {
  # Use custom colors if provided, but fill in missing ones
  st_colors <- custom_st_colors
  # Add auto colors for any missing STs
  missing_st <- setdiff(unique_st, names(custom_st_colors))
  if (length(missing_st) > 0) {
    cat("WARNING: Adding auto colors for missing Sequence Types:", paste(missing_st, collapse=", "), "\n")
    auto_colors <- setNames(
      hcl.colors(length(missing_st), palette = "Pastel 1"),
      missing_st
    )
    st_colors <- c(st_colors, auto_colors)
  }
  cat("Using custom colors for Sequence Type (ST)\n")
} else {
  # Auto-generate colors for all STs
  st_colors <- setNames(
    hcl.colors(length(unique_st), palette = "Pastel 1"),
    unique_st
  )
}

# Print color assignments for verification
cat("\n=== Color Assignments ===\n")
cat("Resistance Profiles:", paste(names(profile_colors), collapse=", "), "\n")
cat("Farms:", paste(names(farm_colors), collapse=", "), "\n")
cat("Diet Stages:", paste(names(diet_colors), collapse=", "), "\n")
cat("Sequence Types:", paste(names(st_colors), collapse=", "), "\n\n")

# ==============================================================================
# 5. CREATE ROW ANNOTATIONS (METADATA FOR ISOLATES)
# ==============================================================================

# Set the order of Diet_Stage levels (Pre-WM, WM, Post-WM)
# This ensures the heatmap rows are split in the correct order
# If your diet stages have different names, update this line:
metadata$Diet_Stage <- factor(metadata$Diet_Stage, 
                              levels = c("Pre-WM", "WM", "Post-WM"))
# Example for different stages:
# levels = c("Baseline", "Treatment", "Recovery")
# levels = c("Day0", "Day7", "Day14", "Day21")

row_ha <- rowAnnotation(
  `Resistance Profile` = metadata$Resistance_Profile,
  Farm = metadata$Farm,
  `Diet Stage` = metadata$Diet_Stage,
  ST = metadata$ST,
  col = list(
    `Resistance Profile` = profile_colors,
    Farm = farm_colors,
    `Diet Stage` = diet_colors,
    ST = st_colors
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.5, "cm"),
  gap = unit(0.5, "mm")
)

# ==============================================================================
# 6. CREATE COLUMN ANNOTATIONS (ANTIMICROBIAL CLASSES)
# ==============================================================================

# Create color scheme for antimicrobial classes
unique_classes <- unique(gene_info$Classification)
class_colors <- setNames(
  rainbow(length(unique_classes)),
  unique_classes
)

# ==============================================================================
# 7. CREATE HEATMAPS - ONE PER FARM
# ==============================================================================

# Define color scheme for gene presence/absence
heatmap_colors <- colorRamp2(c(0, 1), c("white", "navy"))

# Get unique farms
farms <- sort(unique(metadata$Farm))

cat("\n=== Creating heatmaps for", length(farms), "farms ===\n")

# Create a separate heatmap for each farm
for (farm in farms) {
  
  cat("\nProcessing Farm", farm, "...\n")
  
  # Subset data for this farm
  farm_idx <- which(metadata$Farm == farm)
  farm_metadata <- metadata[farm_idx, , drop = FALSE]
  farm_gene_matrix <- gene_matrix[farm_idx, , drop = FALSE]
  
  # Check if farm has any data
  if (nrow(farm_gene_matrix) == 0) {
    cat("  WARNING: No isolates found for Farm", farm, ". Skipping.\n")
    next
  }
  
  cat("  Found", nrow(farm_gene_matrix), "isolates for Farm", farm, "\n")
  
  # FILTER OUT GENES NOT PRESENT IN THIS FARM
  # Keep only genes that are present in at least one isolate
  genes_present <- colSums(farm_gene_matrix) > 0
  farm_gene_matrix <- farm_gene_matrix[, genes_present, drop = FALSE]
  farm_gene_info <- gene_info[genes_present, ]
  
  cat("  Genes present in Farm", farm, ":", sum(genes_present), "out of", length(genes_present), "\n")
  
  # Use the filtered matrix (NOT transposed - isolates on rows, genes on columns)
  
  # Set the order of Diet_Stage levels for this farm
  farm_metadata$Diet_Stage <- factor(farm_metadata$Diet_Stage, 
                                     levels = c("Pre-WM", "WM", "Post-WM"))
  
  # Create ROW annotation for isolates (RP, DS, ST on the left - horizontal bars)
  farm_row_ha <- rowAnnotation(
    RP = farm_metadata$Resistance_Profile,
    DS = farm_metadata$Diet_Stage,
    ST = farm_metadata$ST,
    col = list(
      RP = profile_colors,
      DS = diet_colors,
      ST = st_colors
    ),
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(0.5, "cm"),
    gap = unit(0.5, "mm")
  )
  
  # Create COLUMN annotation for genes (antimicrobial class on top)
  farm_col_ha <- HeatmapAnnotation(
    `Antimicrobial Class` = farm_gene_info$Classification,
    col = list(`Antimicrobial Class` = class_colors),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(0.5, "cm")
  )
  
  # Create the heatmap for this farm (isolates on rows, genes on columns)
  ht <- Heatmap(
    farm_gene_matrix,  # NOT transposed
    name = "Gene Presence",
    
    # Colors
    col = heatmap_colors,
    
    # Row settings (ISOLATES - on the left/Y-axis)
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    row_title = paste("Farm", farm, "- Isolates"),
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_dend_width = unit(2, "cm"),
    row_split = farm_metadata$Diet_Stage,  # Split rows by Diet Stage
    # Alternative row splitting options (uncomment ONE of these to use instead):
    # row_split = farm_metadata$Resistance_Profile,  # Split by Resistance Profile
    # row_split = farm_metadata$ST,  # Split by Sequence Type
    row_gap = unit(0.5, "mm"),
    
    # Column settings (GENES - on top/X-axis)
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    column_title = "Antimicrobial Resistance Genes",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_columns = FALSE,  # Keep genes ordered by class
    show_column_dend = FALSE,
    column_split = farm_gene_info$Classification,  # Split by antimicrobial class
    column_gap = unit(0.5, "mm"),
    
    # Annotations
    left_annotation = farm_row_ha,  # RP, DS, ST on the left
    top_annotation = farm_col_ha,   # Antimicrobial class on top
    
    # Heatmap settings
    border = TRUE,
    rect_gp = gpar(col = "gray90", lwd = 0.5),
    
    # Legend
    heatmap_legend_param = list(
      title = "Gene Presence",
      at = c(0, 1),
      labels = c("Absent", "Present"),
      legend_direction = "horizontal",
      legend_width = unit(4, "cm")
    )
  )
  
  # ==============================================================================
  # 8. DISPLAY AND/OR SAVE THE HEATMAP FOR THIS FARM
  # ==============================================================================
  
  # Display the plot interactively (always shown if running in R/RStudio)
  cat("\n  Displaying plot for Farm", farm, "...\n")
  draw(ht, 
       heatmap_legend_side = "bottom",
       annotation_legend_side = "right",
       merge_legend = TRUE)
  
  # Save files if save_files is TRUE
  if (save_files) {
    # Create filenames with farm number
    pdf_filename <- paste0("antimicrobial_resistance_heatmap_Farm_ncbi_", farm, ".pdf")
    png_filename <- paste0("antimicrobial_resistance_heatmap_Farm_ncbi_", farm, ".png")
    
    # Save as PDF
    pdf(pdf_filename, width = 16, height = 12)
    draw(ht, 
         heatmap_legend_side = "bottom",
         annotation_legend_side = "right",
         merge_legend = TRUE)
    dev.off()
    
    # Save as PNG
    png(png_filename, width = 1600, height = 1200, res = 100)
    draw(ht, 
         heatmap_legend_side = "bottom",
         annotation_legend_side = "right",
         merge_legend = TRUE)
    dev.off()
    
    cat("  - Created:", pdf_filename, "\n")
    cat("  - Created:", png_filename, "\n")
  } else {
    cat("  - Preview mode: Files NOT saved (set save_files <- TRUE to save)\n")
  }
}

cat("\n=== All heatmaps generated successfully! ===\n")

if (save_files) {
  cat("\nOutput files created:\n")
  for (farm in farms) {
    cat("  Farm", farm, ":\n")
    cat("    - antimicrobial_resistance_heatmap_Farm_ncbi_", farm, ".pdf\n", sep="")
    cat("    - antimicrobial_resistance_heatmap_Farm_ncbi_", farm, ".png\n", sep="")
  }
} else {
  cat("\nPREVIEW MODE: Plots displayed but NOT saved.\n")
  cat("To save files, set: save_files <- TRUE\n")
}


# ==============================================================================
# 9. CREATE COMBINED PLOT WITH ALL FARMS (SPLIT BY FARM)
# ==============================================================================

cat("\n=== Creating combined plot with all farms ===\n")

# Filter genes present in at least one farm
genes_present_all <- colSums(gene_matrix) > 0
combined_gene_matrix <- gene_matrix[, genes_present_all, drop = FALSE]
combined_gene_info <- gene_info[genes_present_all, ]

cat("Genes present across all farms:", sum(genes_present_all), "out of", length(genes_present_all), "\n")

# Set the order of Diet_Stage levels
metadata$Diet_Stage <- factor(metadata$Diet_Stage, 
                              levels = c("Pre-WM", "WM", "Post-WM"))

# Convert Farm to factor with explicit ordering (Farm 1, Farm 2, Farm 3, etc.)
metadata$Farm <- factor(metadata$Farm, levels = sort(unique(metadata$Farm)))

# Create row annotation for all isolates (including Farm as a color-coded column)
combined_row_ha <- rowAnnotation(
  Farm = metadata$Farm,  # ADD FARM COLUMN
  RP = metadata$Resistance_Profile,
  DS = metadata$Diet_Stage,
  ST = metadata$ST,
  col = list(
    Farm = farm_colors,  # Color code farms
    RP = profile_colors,
    DS = diet_colors,
    ST = st_colors
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.5, "cm"),
  gap = unit(2, "mm")
)

# Create column annotation for genes
combined_col_ha <- HeatmapAnnotation(
  `Antimicrobial Class` = combined_gene_info$Classification,
  col = list(`Antimicrobial Class` = class_colors),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.5, "cm")
)

# Create the combined heatmap
ht_combined <- Heatmap(
  combined_gene_matrix,
  name = "Gene Presence",
  
  # Colors
  col = heatmap_colors,
  
  # Row settings (ISOLATES)
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  row_title = "Isolates (All Farms)",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows = TRUE,  # Cluster within each farm
  cluster_row_slices = FALSE,  # DON'T reorder the farm groups themselves
  show_row_dend = TRUE,
  row_dend_width = unit(2, "cm"),
  row_split = metadata$Farm,  # SPLIT BY FARM (maintains factor order: 1, 2, 3)
  row_gap = unit(3, "mm"),  # Larger gap between farms
  
  # Column settings (GENES)
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  column_title = "Antimicrobial Resistance Genes",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_split = combined_gene_info$Classification,
  column_gap = unit(2, "mm"),
  
  # Annotations
  left_annotation = combined_row_ha,
  top_annotation = combined_col_ha,
  
  # Heatmap settings
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.5),
  
  # Legend
  heatmap_legend_param = list(
    title = "Gene Presence",
    at = c(0, 1),
    labels = c("Absent", "Present"),
    legend_direction = "horizontal",
    legend_width = unit(4, "cm")
  )
)

# Display the combined plot
cat("\nDisplaying combined plot...\n")
draw(ht_combined, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "right",
     merge_legend = TRUE)

# Save combined plot if save_files is TRUE
if (save_files) {
  # Save as PDF
  pdf("antimicrobial_resistance_heatmap_ALL_FARMSncbi_.pdf", width = 18, height = 14)
  draw(ht_combined, 
       heatmap_legend_side = "bottom",
       annotation_legend_side = "right",
       merge_legend = TRUE)
  dev.off()
  
  # Save as PNG
  png("antimicrobial_resistance_heatmap_ALL_FARMSncbi_.png", width = 1800, height = 1400, res = 100)
  draw(ht_combined, 
       heatmap_legend_side = "bottom",
       annotation_legend_side = "right",
       merge_legend = TRUE)
  dev.off()
  
  cat("\nCombined plot saved:\n")
  cat("  - antimicrobial_resistance_heatmap_ALL_FARMSncbi_.pdf\n")
  cat("  - antimicrobial_resistance_heatmap_ALL_FARMSncbi_.png\n")
} else {
  cat("\nPREVIEW MODE: Combined plot displayed but NOT saved.\n")
}


cat("\n=== ALL DONE! ===\n")

