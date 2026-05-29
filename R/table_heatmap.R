library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

source("utils.R")


# parameter to change manually
target_param <- "intensity" 

base_dirs <- c(
  "D:/240417_W1_FPs",
  "D:/240924_W1_FPs",
  "D:/241113_W1_FPs"
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

raw_data_all <- data.frame()

cat("Beginning of data extraction...\n")

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
          extracted_values <- (df_poca$intensity / df_poca$`total ON`) * (3.6 / 300)
        } else if (target_param == "intensity") {
          extracted_values <- df_poca$intensity * (3.6 / 300)
        } else {
          extracted_values <- df_poca[[target_param]]
        }
        
        temp_data <- data.frame(
          Protein = prot_name,
          Plate = basename(plate_dir),
          Well = well,
          value = extracted_values
        )
        raw_data_all <- bind_rows(raw_data_all, temp_data)
      }
    }
  }
}

cat("Extraction finished. Computing matrix...\n")

# reference is meos3.2
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

df_heatmap <- melt(log2_matrix, na.rm = TRUE)
colnames(df_heatmap) <- c("Protein_Y", "Protein_X", "Log2_Ratio")

df_heatmap$Protein_X <- factor(df_heatmap$Protein_X, levels = names(proteins_map))
df_heatmap$Protein_Y <- factor(df_heatmap$Protein_Y, levels = names(proteins_map))


# scale calculation
max_abs_log2 <- max(abs(df_heatmap$Log2_Ratio), na.rm = TRUE)
dynamic_limit <- max_abs_log2 * 1.05
ratio_max <- round(2^dynamic_limit, 2)
ratio_min <- round(2^(-dynamic_limit), 2)


p_heatmap <- ggplot(df_heatmap, aes(x = Protein_X, y = Protein_Y, fill = Log2_Ratio)) +
  geom_tile(color = "white", linewidth = 0.5) + 
  
  scale_fill_distiller(
    palette = "RdBu", 
    direction = -1, 
    limits = c(-dynamic_limit, dynamic_limit),
    breaks = c(-dynamic_limit, 0, dynamic_limit),
    labels = c(paste0("x ", ratio_min), "x 1", paste0("x ", ratio_max)),
    na.value = "whitesmoke"
  ) +
  
  scale_y_discrete(limits = rev) +
  
  labs(title = paste("Protein Comparison Matrix:", target_param),
       fill = "Ratio") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "whitesmoke", color = NA)
  ) +
  coord_fixed()

print(p_heatmap)


ggsave(filename = paste0("comparison_matrix_", target_param, ".png"), 
       plot = p_heatmap, 
       width = 8, 
       height = 6, 
       dpi = 300, 
       bg = "white")