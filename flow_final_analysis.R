library(ggplot2)

#####################################################################
# Reusable CONSORT-style attrition flow chart
#
# `stages` is a data.frame with columns:
#   label     - text describing the stage
#   n_isolates
#   n_samples
#####################################################################

draw_flowchart <- function(stages, title = "") {
  
  n_stages <- nrow(stages)
  box_height <- 1
  gap <- 1.4   # vertical spacing between boxes
  box_width <- 8
  
  stages$y <- rev(seq_len(n_stages)) * gap
  stages$x <- 0
  
  stages$text <- sprintf("%s\n%d isolates | %d samples",
                         stages$label, stages$n_isolates, stages$n_samples)
  
  # Loss annotations (difference between consecutive stages)
  loss_labels <- data.frame(
    x = 4.3,
    y = (stages$y[-n_stages] + stages$y[-1]) / 2,
    label = sprintf("- %d isolates\n- %d samples\nexcluded",
                    stages$n_isolates[-n_stages] - stages$n_isolates[-1],
                    stages$n_samples[-n_stages] - stages$n_samples[-1])
  )
  
  p <- ggplot() +
    # Boxes
    geom_rect(data = stages,
              aes(xmin = x - box_width/2, xmax = x + box_width/2,
                  ymin = y - box_height/2, ymax = y + box_height/2),
              fill = "white", color = "black") +
    geom_text(data = stages, aes(x = x, y = y, label = text), size = 3.5, lineheight = 0.9) +
    # Arrows between boxes
    geom_segment(data = stages[-n_stages, ],
                 aes(x = 0, xend = 0, y = y - box_height/2, yend = y - gap + box_height/2),
                 arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5) +
    # Loss annotations to the side of each arrow
    geom_text(data = loss_labels, aes(x = x, y = y, label = label),
              size = 3, color = "firebrick", hjust = 0, lineheight = 0.9) +
    coord_cartesian(clip = "off") +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.margin = margin(10, 60, 10, 10))
  
  p
}

#####################################################################
# Example usage - E. coli pipeline
# (uses the objects already created in ecoli_full_pipeline_valid_samples.R:
#  esch_all, esch, amr)
#####################################################################

ecoli_stages <- data.frame(
  label = c(
    "All E. coli isolates\n(NGIsolates.xlsx)",
    "After filtering to final/faecal\nsamples (calves_per_farm_final.xlsx)",
    "Final resistance analysis set"
  ),
  n_isolates = c(nrow(esch_all), nrow(esch), nrow(amr)),
  n_samples  = c(n_distinct(esch_all$FKSampleID), n_distinct(esch$FKSampleID), n_distinct(amr$FKSampleID))
)

p_flow_ecoli <- draw_flowchart(ecoli_stages, title = "E. coli sample/isolate attrition")
print(p_flow_ecoli)
ggsave("ecoli_attrition_flowchart.png", p_flow_ecoli, width = 8, height = 6, dpi = 300)

#####################################################################
# Example usage - Enterococcus pipeline
# (uses the objects from enterococcus_full_pipeline_valid_samples.R:
#  entero_all, entero, amr)
#####################################################################

 entero_stages <- data.frame(
   label = c(
     "All Enterococcus (faecalis/faecium)\nisolates (NGIsolates.xlsx)",
     "After filtering to final/faecal\nsamples (calves_per_farm_final.xlsx)",
          "Final resistance analysis set"
   ),
   n_isolates = c(nrow(entero_all), nrow(entero), nrow(amr)),
   n_samples  = c(n_distinct(entero_all$FKSampleID), n_distinct(entero$FKSampleID), n_distinct(amr$FKSampleID))
 )

 p_flow_entero <- draw_flowchart(entero_stages, title = "Enterococcus sample/isolate attrition")
 print(p_flow_entero)
 ggsave("enterococcus_attrition_flowchart.png", p_flow_entero, width = 8, height = 6, dpi = 300)