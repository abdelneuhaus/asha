from asha.setup_plate import setup_plate
from asha.src.io_utils import get_PALMTracer_files, get_poca_files
from plot_gradient_spycatcher import plot_gradient_spycatcher
from plot_gradient_meos import plot_gradient_meos
from asha.src.io_utils import remove_files

plate_gamme1 = "D:/ANALYSIS_PAPER/test_gammes/MEOS"   # 250213_HCS_gamme, 250205_HCS_gamme, 250123_HCS_gamme, 250130_HCS_gamme
display_mode = "random" # random or both
list_of_pt_files = get_PALMTracer_files(plate_gamme1)
list_of_poca_files = get_poca_files(plate_gamme1)
# plot_gradient_spycatcher(list_of_pt_files, list_of_poca_files, display_mode)
plot_gradient_meos(list_of_pt_files, list_of_poca_files, display_mode)

# remove_files("D:/ANALYSIS_PAPER/250123_HCS_gamme")
# remove_files("D:/ANALYSIS_PAPER/250130_HCS_gamme")
# remove_files("D:/ANALYSIS_PAPER/250205_HCS_gamme")
# remove_files("D:/ANALYSIS_PAPER/test_gammes/MEOS_BIS")
