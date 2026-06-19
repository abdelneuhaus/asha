library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)

source("utils.R")


base_dir <- "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs"
params_to_plot <- c("photon_loc", "intensity", "blinks", "avg_on", "avg_off", "total ON") 

protein_db <- list(
  "mEos4b-L93M" = c('C3', 'D9', 'E3'),
  "mEos3.2"     = c('C4', 'D8', 'E4'),
  "mEosEM"      = c('C5', 'D7', 'E5'),
  "pcSTAR"      = c('C6', 'D6', 'E6'),
  "Dendra2"     = c('C7', 'D5', 'E7'),
  "Maple3"      = c('C8', 'D4', 'E8'),
  "mEos4b_"     = c('C9', 'D3', 'E9')
)

well_colors <- c("Well 1" = "#6BAED6", "Well 2" = "#3182BD", "Well 3" = "#08519C")

custom_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none"
  )

for (param in params_to_plot) {
  
  plot_list <- list()
  message(paste("Génération de la figure pour :", param))
  
  for (prot_name in names(protein_db)) {
    target_wells <- protein_db[[prot_name]]
    mean_data_prot <- data.frame()
    
    for (w_idx in seq_along(target_wells)) {
      well <- target_wells[w_idx]
      well_dir <- file.path(base_dir, well)
      poca_files <- get_poca_files(well_dir)
      
      if (length(poca_files) == 0) next
      
      for (f in seq_along(poca_files)) {
        df_poca <- read_poca_files(poca_files[f])
        
        if (param == "photon_loc") {
          if (all(c("intensity", "total ON") %in% colnames(df_poca))) {
            extracted_values <- (df_poca$intensity / df_poca$`total ON`) * (3.6 / 300)
          } else next
          
        } else if (param == "avg_on") {
          if (all(c("total ON", "# seq ON") %in% colnames(df_poca))) {
            extracted_values <- ifelse(df_poca$`# seq ON` != 0, df_poca$`total ON` / df_poca$`# seq ON`, 0)
          } else next
          
        } else if (param == "avg_off") {
          if (all(c("total OFF", "# seq OFF") %in% colnames(df_poca))) {
            extracted_values <- ifelse(df_poca$`# seq OFF` != 0, df_poca$`total OFF` / df_poca$`# seq OFF`, 0)
          } else next
          
        } else {
          if (param %in% colnames(df_poca)) {
            extracted_values <- df_poca[[param]]
          } else next
        }
        
        mean_data_prot <- bind_rows(mean_data_prot, data.frame(
          Well_Original = well,
          Well_Label = paste("Well", w_idx),
          mean_value = mean(extracted_values, na.rm = TRUE)
        ))
      }
    }
    
    if (nrow(mean_data_prot) == 0) {
      warning(paste("No data for", prot_name, "(parameter:", param, ")"))
      next
    }
    
    mean_data_prot$Well_Label <- factor(mean_data_prot$Well_Label, levels = c("Well 1", "Well 2", "Well 3"))
    
    p <- ggplot(mean_data_prot, aes(x = Well_Label, y = mean_value, fill = Well_Label)) +
      stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7, alpha = 0.8) +
      stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1) +
      geom_jitter(width = 0.05, shape = 21, fill = "white", color = "black", size = 2, alpha = 0.9) +
      scale_fill_manual(values = well_colors) +
      stat_compare_means(method = "kruskal.test", label.y = max(mean_data_prot$mean_value, na.rm = TRUE) * 1.15) + 
      labs(title = prot_name, y = param, x = "") +
      custom_theme
    
    plot_list[[prot_name]] <- p
  }
  
  
  if (length(plot_list) > 0) {
    final_plot <- wrap_plots(plot_list, ncol = 4, nrow = 2) +
      plot_annotation(
        title = paste("Photophysical Parameter :", toupper(param)),
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
    
    print(final_plot)
    
    filename <- paste0("figure_wells_", param, ".png")
    ggsave(filename, plot = final_plot, width = 16, height = 8, dpi = 300)
    message("saved : ", filename, "\n")
  }
}