library(dplyr)
library(ggpubr)
library(ggplot2)
library(patchwork)

source("utils.R")

# parameter to change manually
target_param <- "photon_loc" 

# title to change manually
y_axis_label <- "Number of photons\nper localization" 
y_limits <- c(45, 210) 

base_dirs <- c(
  "D:/240417_W1_FPs",
  "D:/240924_W1_FPs",
  "D:/241113_W1_FPs"
)
plate_names <- c("Plate 1", "Plate 2", "Plate 3")
target_wells <- c("E3", "D9", "C3")

raw_data_fov <- data.frame() 
mean_data_all <- data.frame()

for (p in seq_along(base_dirs)) {
  plate_dir <- base_dirs[p]
  plate_name <- plate_names[p]
  
  for (well in target_wells) {
    well_dir <- file.path(plate_dir, well)
    
    poca_files <- get_poca_files(well_dir)
    
    if (length(poca_files) == 0) next
    
    for (f in seq_along(poca_files)) {
      df_poca <- read_poca_files(poca_files[f])
      
      if (target_param == "photon_loc") {
        if ("intensity" %in% colnames(df_poca) && "total ON" %in% colnames(df_poca)) {
          extracted_values <- (df_poca$intensity / df_poca$`total ON`) * (3.6 / 300)
        } else {
          warning(paste("'intensity' or 'total ON' columns missing in:", poca_files[f]))
          next
        }
      } else {
        if (target_param %in% colnames(df_poca)) {
          extracted_values <- df_poca[[target_param]]
        } else {
          warning(paste("Missing column (", target_param, ") in file:", poca_files[f]))
          next
        }
      }
      
      if (p == 1 && well == target_wells[1]) {
        temp_raw <- data.frame(
          FOV = paste("FOV", f),
          value = extracted_values
        )
        raw_data_fov <- bind_rows(raw_data_fov, temp_raw)
      }
      
      temp_mean <- data.frame(
        Plate = plate_name,
        Well = well,
        FOV = paste("FOV", f),
        mean_value = mean(extracted_values, na.rm = TRUE)
      )
      mean_data_all <- bind_rows(mean_data_all, temp_mean)
    }
  }
}


mean_data_all$Well <- factor(mean_data_all$Well, levels = target_wells)
levels(mean_data_all$Well) <- c("Well 1", "Well 2", "Well 3")

color_well1 <- "#6BAED6"
color_well2 <- "#3182BD"
color_well3 <- "#08519C"
well_colors <- c("Well 1" = color_well1, "Well 2" = color_well2, "Well 3" = color_well3)



custom_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

# fov replicates 
p1 <- ggplot(raw_data_fov, aes(x = FOV, y = value)) +
  geom_boxplot(fill = color_well1, color = "black", outlier.shape = NA, width = 0.5, alpha = 0.7) +
  labs(title = "Field-of-view Replicates", 
       y = y_axis_label, 
       x = "") +
  custom_theme

if (!is.null(y_limits)) { p1 <- p1 + coord_cartesian(ylim = y_limits) }

# well replicates 
data_p2 <- mean_data_all %>% filter(Plate == "Plate 1")
p2 <- ggplot(data_p2, aes(x = Well, y = mean_value, fill = Well)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7, alpha = 0.8) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.05, shape = 21, fill = "white", color = "black", size = 2, alpha = 0.9) +
  scale_fill_manual(values = well_colors) +
  stat_compare_means(method = "kruskal.test", label.y = max(data_p2$mean_value, na.rm = TRUE) * 1.1) + 
  labs(title = "Well Replicates", y = "", x = "") +
  custom_theme

# plate replicates 
p3 <- ggplot(mean_data_all, aes(x = Plate, y = mean_value, fill = Plate)) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.5, alpha = 0.5) +
  geom_jitter(width = 0.1, shape = 21, fill = "white", color = "black", size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c("#377eb8", "#ff7f00", "#4daf4a")) + 
  labs(title = "Plate Replicates", y = "", x = "") +
  custom_theme

if (!is.null(y_limits)) { p3 <- p3 + coord_cartesian(ylim = y_limits) }

final_plot <- p1 + p2 + p3
print(final_plot)
ggsave("figure_replicates_carrés.png", plot = final_plot, width = 12, height = 4.5, dpi = 300)
