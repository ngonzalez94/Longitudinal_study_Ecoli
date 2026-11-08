# ==============================================================================
# Plasmid-ARG Co-localisation Analysis
# Matches PlasmidFinder contigs with NCBI ARG contigs to identify
# resistance genes likely carried on plasmids
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ==============================================================================
# 1. FILE PATHS — UPDATE THESE
# ==============================================================================

plasmidfinder_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/PLASMIDFINDER.xlsx"

ncbi_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/results_ncbi.xlsx"
ncbi_sheet <- 1   # sheet with NCBI results (long format with SEQUENCE column)
# NOTE: this needs to be the RAW abricate output (with SEQUENCE column),
# not the summarised wide matrix. Update the path/sheet if stored separately.

save_files <- TRUE

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

pf  <- read_excel(plasmidfinder_path)
arg <- read_excel(ncbi_path, sheet = ncbi_sheet)

# Clean isolate names
pf$`#FILE`  <- gsub("_contig\\.fa$", "", as.character(pf$`#FILE`))
arg$`#FILE` <- gsub("_contig\\.fa$", "", as.character(arg$`#FILE`))

cat("=== Data loaded ===\n")
cat("PlasmidFinder rows :", nrow(pf),  "\n")
cat("NCBI ARG rows      :", nrow(arg), "\n")
cat("Unique isolates (PF) :", length(unique(pf$`#FILE`)),  "\n")
cat("Unique isolates (ARG):", length(unique(arg$`#FILE`)), "\n\n")

# ==============================================================================
# 3. FIND CONTIG OVERLAPS
#    A gene is "plasmid-associated" if its SEQUENCE (contig) matches a contig
#    that also carries a plasmid replicon in the same isolate
# ==============================================================================

# Build a lookup: isolate + contig -> plasmid(s) found on that contig
plasmid_contigs <- pf %>%
  select(Isolate  = `#FILE`,
         Contig   = SEQUENCE,
         Plasmid  = GENE,
         Plasmid_pct_id = `%IDENTITY`) %>%
  distinct()

# Join ARG data with plasmid contig lookup on isolate + contig
arg_clean <- arg %>%
  select(Isolate  = `#FILE`,
         Contig   = SEQUENCE,
         ARG      = GENE,
         ARG_pct_id = `%IDENTITY`,
         Product  = PRODUCT,
         Resistance_class = RESISTANCE) %>%
  distinct()

plasmid_args <- arg_clean %>%
  inner_join(plasmid_contigs, by = c("Isolate", "Contig"),
             relationship = "many-to-many") %>%
  arrange(Isolate, Contig, ARG)

cat("=== Overlap Results ===\n")
cat("ARG-plasmid contig matches found:", nrow(plasmid_args), "\n")
cat("Unique isolates with plasmid-ARGs:", length(unique(plasmid_args$Isolate)), "\n")
cat("Unique ARGs on plasmid contigs   :", length(unique(plasmid_args$ARG)), "\n")
cat("Unique plasmids carrying ARGs    :", length(unique(plasmid_args$Plasmid)), "\n\n")

# Print summary table
cat("=== ARGs found on plasmid-carrying contigs ===\n")
summary_tbl <- plasmid_args %>%
  group_by(ARG, Plasmid, Resistance_class) %>%
  summarise(n_isolates = n_distinct(Isolate),
            Isolates   = paste(sort(unique(Isolate)), collapse = ", "),
            .groups    = "drop") %>%
  arrange(desc(n_isolates), ARG)
print(summary_tbl, n = Inf)

# ==============================================================================
# 4. SAVE FULL RESULTS TABLE
# ==============================================================================

if (save_files) {
  write.csv(plasmid_args,
            "plasmid_ARG_colocalization_full.csv",
            row.names = FALSE)
  write.csv(summary_tbl,
            "plasmid_ARG_colocalization_summary.csv",
            row.names = FALSE)
  cat("\nSaved:\n")
  cat("  plasmid_ARG_colocalization_full.csv\n")
  cat("  plasmid_ARG_colocalization_summary.csv\n")
}

# ==============================================================================
# 5. HEATMAP: isolates (rows) vs ARGs (columns)
#    Cell colour = which plasmid replicon carries that ARG in that isolate
#    White = ARG not on a plasmid contig in that isolate
# ==============================================================================

# Build presence matrix: for each isolate x ARG, record plasmid name (or NA)
# If multiple plasmids on the same contig, keep the first for simplicity
heat_data <- plasmid_args %>%
  group_by(Isolate, ARG) %>%
  summarise(Plasmid = first(Plasmid), .groups = "drop")

# All isolates that have ANY plasmid-ARG
all_isolates <- sort(unique(plasmid_args$Isolate))
all_args     <- sort(unique(plasmid_args$ARG))

# Binary presence matrix (1 = ARG on a plasmid contig, 0 = not)
presence_mat <- matrix(0L,
                       nrow = length(all_isolates),
                       ncol = length(all_args),
                       dimnames = list(all_isolates, all_args))

for (i in seq_len(nrow(heat_data))) {
  iso <- heat_data$Isolate[i]
  gen <- heat_data$ARG[i]
  if (iso %in% rownames(presence_mat) && gen %in% colnames(presence_mat)) {
    presence_mat[iso, gen] <- 1L
  }
}

if (nrow(presence_mat) == 0 || ncol(presence_mat) == 0) {
  cat("\nNo plasmid-ARG overlaps found — heatmap skipped.\n")
} else {
  
  # Color: white = absent, coral = present on plasmid contig
  heat_colors <- colorRamp2(c(0, 1), c("white", "#E31A1C"))
  
  ht_plasmid <- Heatmap(
    presence_mat,
    name                = "On plasmid contig",
    col                 = heat_colors,
    show_heatmap_legend = FALSE,
    
    row_names_side  = "left",
    row_names_gp    = gpar(fontsize = 8),
    row_title       = "Isolates",
    row_title_gp    = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows    = TRUE,
    show_row_dend   = TRUE,
    
    column_names_side = "bottom",
    column_names_gp   = gpar(fontsize = 8),
    column_names_rot  = 45,
    column_title      = "Resistance Genes on Plasmid-Carrying Contigs",
    column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
    cluster_columns   = TRUE,
    show_column_dend  = TRUE,
    
    border  = TRUE,
    rect_gp = gpar(col = "gray90", lwd = 0.5),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (presence_mat[i, j] == 1) {
        # Show the plasmid name inside the cell
        iso <- rownames(presence_mat)[i]
        gen <- colnames(presence_mat)[j]
        plas <- heat_data$Plasmid[heat_data$Isolate == iso &
                                    heat_data$ARG    == gen][1]
        if (!is.na(plas)) {
          grid.text(plas, x, y,
                    gp = gpar(fontsize = 5, col = "white", fontface = "bold"))
        }
      }
    }
  )
  
  presence_legend <- Legend(
    title     = "ARG on plasmid contig",
    labels    = c("No", "Yes"),
    legend_gp = gpar(fill = c("white", "#E31A1C"), col = "gray60"),
    border    = "gray60"
  )
  
  draw(ht_plasmid,
       annotation_legend_side = "right",
       annotation_legend_list = list(presence_legend))
  
  if (save_files) {
    pdf("plasmid_ARG_heatmap.pdf", width = 14, height = 10)
    draw(ht_plasmid, annotation_legend_side = "right",
         annotation_legend_list = list(presence_legend))
    dev.off()
    png("plasmid_ARG_heatmap.png", width = 1400, height = 1000, res = 100)
    draw(ht_plasmid, annotation_legend_side = "right",
         annotation_legend_list = list(presence_legend))
    dev.off()
    cat("\nSaved: plasmid_ARG_heatmap.pdf & plasmid_ARG_heatmap.png\n")
  }
}

# ==============================================================================
# 6. BUBBLE PLOT: plasmid type vs ARG, bubble size = number of isolates
# ==============================================================================

if (nrow(plasmid_args) > 0) {
  
  bubble_data <- plasmid_args %>%
    group_by(Plasmid, ARG, Resistance_class) %>%
    summarise(n_isolates = n_distinct(Isolate), .groups = "drop")
  
  p_bubble <- ggplot(bubble_data,
                     aes(x = ARG, y = Plasmid, size = n_isolates,
                         fill = Resistance_class)) +
    geom_point(shape = 21, colour = "gray40", alpha = 0.85) +
    scale_size_continuous(name = "No. isolates", range = c(3, 12)) +
    scale_fill_brewer(name = "Resistance class", palette = "Set1",
                      na.value = "grey70") +
    labs(title = "Plasmid replicons co-localised with resistance genes",
         x = "Resistance gene", y = "Plasmid replicon") +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y  = element_text(size = 8),
          panel.grid.major = element_line(colour = "gray90"),
          plot.title   = element_text(face = "bold", size = 12))
  
  print(p_bubble)
  
  if (save_files) {
    ggsave("plasmid_ARG_bubbleplot.pdf", p_bubble,
           width = max(8, length(unique(bubble_data$ARG)) * 0.5 + 3),
           height = max(6, length(unique(bubble_data$Plasmid)) * 0.4 + 3),
           limitsize = FALSE)
    ggsave("plasmid_ARG_bubbleplot.png", p_bubble,
           width = max(8, length(unique(bubble_data$ARG)) * 0.5 + 3),
           height = max(6, length(unique(bubble_data$Plasmid)) * 0.4 + 3),
           dpi = 150, limitsize = FALSE)
    cat("Saved: plasmid_ARG_bubbleplot.pdf & plasmid_ARG_bubbleplot.png\n")
  }
}

cat("\n=== ALL DONE! ===\n")
