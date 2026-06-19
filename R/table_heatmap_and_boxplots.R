library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

source("utils.R")

# parameters
target_params <- c("intensity", "photon_loc", "total ON", "avg_on", "avg_off", "blinks") 

base_dirs <- c(
  # "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs"
  # "D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs"
  "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs"
)

proteins_map <- list(
  "mEos3.2" = c("C4", "D8", "E4"),
  "mEos4b" = c("C9", "D3", "E9"),
  "mEos4b-L93M" = c("C3", "D9", "E3"),
  "pcSTAR" = c("C6", "D6", "E6"),
  "mEosEM" = c("C5", "D7", "E5"),
  "Dendra2" = c("C7", "D5", "E7"),
  "mMaple3" = c("C8", "D4", "E8")
)

plot_rows_list <- list()

# main loop per protein
for (target_param in target_params) {
  cat(paste("\nAnalyzing parameter :", target_param, "---\n"))
  raw_data_all <- data.frame()
  
  for (plate_dir in base_dirs) {
    for (prot_name in names(proteins_map)) {
      wells <- proteins_map[[prot_name]]
      for (well in wells) {
        well_dir <- file.path(plate_dir, well)
        poca_files <- get_poca_files(well_dir)
        if (length(poca_files) == 0) next
        
        for (f in seq_along(poca_files)) {
          df_poca <- read_poca_files(poca_files[f])
          
          if (target_param == "photon_loc") {
            if (all(c("intensity", "total ON") %in% colnames(df_poca))) {
              extracted_values <- (df_poca$intensity / df_poca$`total ON`) * (3.6 / 300)
            } else next
            
          } else if (target_param == "intensity") {
            if ("intensity" %in% colnames(df_poca)) {
              extracted_values <- df_poca$intensity * (3.6 / 300)
            } else next
            
          } else if (target_param == "avg_on") {
            if (all(c("total ON", "# seq ON") %in% colnames(df_poca))) {
              extracted_values <- ifelse(df_poca$`# seq ON` != 0, df_poca$`total ON` / df_poca$`# seq ON`, 0)
            } else next
            
          } else if (target_param == "avg_off") {
            if (all(c("total OFF", "# seq OFF") %in% colnames(df_poca))) {
              extracted_values <- ifelse(df_poca$`# seq OFF` != 0, df_poca$`total OFF` / df_poca$`# seq OFF`, 0)
            } else next
            
          } else {
            if (target_param %in% colnames(df_poca)) {
              extracted_values <- df_poca[[target_param]]
            } else next
          }
          
          extracted_values <- extracted_values[is.finite(extracted_values)]
          
          temp_data <- data.frame(
            Protein = prot_name,
            Plate = basename(plate_dir),
            value = extracted_values
          )
          raw_data_all <- bind_rows(raw_data_all, temp_data)
        }
      }
    }
  }
  
  if (nrow(raw_data_all) == 0) {
    cat(paste("No data for", target_param, "- skipped\n"))
    next
  }
  
  # matrix norm
  ref_data <- raw_data_all %>%
    filter(Protein == "mEos3.2") %>%
    group_by(Plate) %>%
    summarize(ref_mean = mean(value, na.rm = TRUE))
  
  normalized_data <- raw_data_all %>%
    left_join(ref_data, by = "Plate") %>%
    mutate(norm_value = value / ref_mean) 
  
  summary_data <- normalized_data %>%
    group_by(Protein) %>%
    summarize(mean_val = mean(norm_value, na.rm = TRUE))
  
  protein_means <- summary_data$mean_val
  names(protein_means) <- summary_data$Protein
  protein_means <- protein_means[names(proteins_map)]
  
  ratio_matrix <- outer(protein_means, protein_means, FUN = "/")
  log2_matrix <- log2(ratio_matrix)
  
  df_heatmap <- melt(log2_matrix, na.rm = TRUE)
  colnames(df_heatmap) <- c("Protein_Y", "Protein_X", "Log2_Ratio")
  
  df_heatmap$Protein_X <- factor(df_heatmap$Protein_X, levels = names(proteins_map))
  df_heatmap$Protein_Y <- factor(df_heatmap$Protein_Y, levels = names(proteins_map))
  
  max_abs_log2 <- max(abs(df_heatmap$Log2_Ratio), na.rm = TRUE)
  dynamic_limit <- max_abs_log2 * 1.05
  ratio_max <- round(2^dynamic_limit, 2)
  ratio_min <- round(2^(-dynamic_limit), 2)
  
  p_heatmap <- ggplot(df_heatmap, aes(x = Protein_X, y = Protein_Y, fill = Log2_Ratio)) +
    geom_tile(color = "white", linewidth = 0.5) + 
    scale_fill_distiller(
      palette = "RdBu", direction = -1, 
      limits = c(-dynamic_limit, dynamic_limit), 
      breaks = c(-dynamic_limit, 0, dynamic_limit), 
      labels = c(paste0("x ", ratio_min), "x 1", paste0("x ", ratio_max)), 
      na.value = "whitesmoke"
    ) +
    scale_y_discrete(limits = rev) +
    labs(title = paste("Matrix:", toupper(target_param)), fill = "Ratio") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "whitesmoke", color = NA)
    ) +
    coord_fixed()
  
  # boxplots
  raw_data_all$Protein <- factor(raw_data_all$Protein, levels = names(proteins_map))
  
  y_max_zoom <- quantile(raw_data_all$value, 0.99, na.rm = TRUE)
  
  p_boxplot <- ggplot(raw_data_all, aes(x = Protein, y = value, fill = Protein)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.6) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste("Global Distribution:", toupper(target_param)), y = target_param, x = "") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.title.y = element_text(size = 12),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0, y_max_zoom)) 
  
  p_row <- p_heatmap + p_boxplot + plot_layout(widths = c(1, 1.5))
  plot_rows_list[[target_param]] <- p_row
}

cat("\nAssemblage final de toutes les lignes...\n")


final_mega_plot <- wrap_plots(plot_rows_list, ncol = 1)
total_height <- length(plot_rows_list) * 5
ggsave(filename = "mega_figure.png", 
       plot = final_mega_plot, 
       width = 16, 
       height = total_height, 
       dpi = 300, 
       bg = "white",
       limitsize = FALSE)

cat("Finished!'\n")