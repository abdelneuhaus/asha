from setup_plate import setup_plate
from utils import remove_files
from boxplots_inside_a_well import boxplots_inside_a_well
from boxplots_inside_the_plate import boxplots_inside_the_plate
from boxplots_between_plate import boxplots_between_plate, generate_data_plate_boxplot
from plot_heatmap_one_parameter import plot_heatmap_one_parameter


# plate0 = "D:/ANALYSIS_PAPER/240313_W1_FPs"
# plate1 = "D:/ANALYSIS_PAPER/240326_W1_FPs"
# plate2 = "D:/ANALYSIS_PAPER/240417_W1_FPs"
# plate3 = "D:/ANALYSIS_PAPER/240924_W1_FPs"
# plate4 = "D:/ANALYSIS_PAPER/241113_W1_FPs"
# plate5 = "D:/ANALYSIS_PAPER/250326_W1_FPs"

plate2bis = "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs"
plate3bis = "D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs"
plate4bis = "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs"

prefix = 'SR_001.MIA'
group_poca_files, group_palm_tracer_files = setup_plate(plate3bis, prefix)

protein = "mEos3.2"
parameter = "Integrated_Intensity"  # or "intensity", "avg_on", "photon_loc", "total ON", "blinks", "Integrated_Intensity"

# List of target proteins to analyze: meos4b-L93M, mEos3.2, mEosEM, pcSTAR, Dendra2, mMaple3, mEos4b
# boxplots_inside_a_well(protein, group_poca_files, group_palm_tracer_files)


# Boxplots inside a plate for a specific parameter inside the same plate
# boxplots_inside_the_plate(group_poca_files, group_palm_tracer_files)

# Boxplots inside a plate for a specific protein and parameter, and comparison between plates
# boxplot_plate0 = generate_data_plate_boxplot(plate0, protein, parameter, prefix)
# boxplot_plate1 = generate_data_plate_boxplot(plate1, protein, parameter, prefix)
# boxplot_plate2 = generate_data_plate_boxplot(plate2, protein, parameter, prefix)
# boxplot_plate3 = generate_data_plate_boxplot(plate3, protein, parameter, prefix)
# boxplot_plate4 = generate_data_plate_boxplot(plate4, protein, parameter, prefix)
# # boxplot_plate5 = generate_data_plate_boxplot(plate5, protein, parameter, prefix)
# boxplot_plate2bis = generate_data_plate_boxplot(plate2bis, protein, parameter, prefix)
# boxplot_plate3bis = generate_data_plate_boxplot(plate3bis, protein, parameter, prefix)
# boxplot_plate4bis = generate_data_plate_boxplot(plate4bis, protein, parameter, prefix)
# # # boxplots_between_plate([boxplot_plate0, boxplot_plate1, boxplot_plate2, boxplot_plate3, boxplot_plate4, boxplot_plate5])
# boxplots_between_plate([boxplot_plate2bis, boxplot_plate3bis, boxplot_plate4bis])


# # Plate heatmap for a specific parameter
# plot_heatmap_one_parameter(plate2bis, parameter, prefix)
# plot_heatmap_one_parameter(plate3bis, parameter, prefix)
# plot_heatmap_one_parameter(plate4bis, parameter, prefix)


# Remove files if needed
# remove_files(plate0)
# remove_files(plate1)
# remove_files(plate2)
# remove_files(plate3)
# remove_files(plate4)
# remove_files(plate5)
# remove_files("D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs")
# remove_files("D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs")
# remove_files("D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs")
# remove_files("D:/ANALYSIS_PAPER/new_threshold/240313_W1_FPs")
# remove_files("D:/ANALYSIS_PAPER/new_threshold/240326_W1_FPs")
# remove_files("D:/ANALYSIS_PAPER/new_threshold/250326_W1_FPs")
