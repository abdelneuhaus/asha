from tkinter import *					
from tkinter import ttk
from tkinter.filedialog import askdirectory

from PAPER.utils import read_poca_files, get_poca_files, get_PALMTracer_files, get_csv_poca_frame_files, get_csv_poca_intensity_files, get_csv_poca_sigma_files
from do_analysis_for_one_acquisition import do_analysis_for_one_acquisition
from do_cumulative_number_clusters import do_cumulative_number_clusters
from asha.heatmap_single_parameters import do_96heatmap_one_photophysics_parameter
from time_distribution_same_well import time_distribution_same_well
from time_distribution_same_row import time_distribution_same_row
from asha.preprocessing import fusion
from fov_photophysics_distribution import fov_photophysics_distribution
from do_fps_heatmap import do_fps_heatmap, generate_fps_heatmap

import os
import statistics
import pandas as pd
import numpy as np


class MyWindow:
    
    def __init__(self, root):
        window = ttk.Notebook(root)

        # Tab 1 : Load files, etc...
        tab1 = ttk.Frame(window)
        window.add(tab1, text ='Files & Informations')
        # Tab 2 : Global Experiment Stats (Density, Photophysics Plot)
        tab2 = ttk.Frame(window)
        window.add(tab2, text ='Global Stats Analysis')
        # Tab 3 : 8 well Plate Analysis
        tab3 = ttk.Frame(window)
        window.add(tab3, text ='8-Well Plate Analysis')
        # Tab 4 : 96 well Plate Analysis (HCS)
        tab4 = ttk.Frame(window)
        window.add(tab4, text ='96-Well Plate Analysis')
        window.pack(expand = 1, fill ="both")

        # ------- INFORMATIONS TAB -------
        self.load_data_bool = BooleanVar()
        self.load_data_bool.set(False)
        self.load_data = Button(tab1, text='Select Experiment', command=self.load_molecule_data)
        self.load_data.grid(row=0, column=0, columnspan=7, sticky="WE", pady=3, ipadx=1, padx=5)
        self.repertory_path = None
        self.index = None
        self.files = None
        self.laser = ['561.PT', '561-405.PT']
        
        self.isPT_bool = BooleanVar()
        self.isPT_bool.set(True)
        self.isPT = Checkbutton(tab1, text='PALMTracer experiment (After Loading)', variable=self.isPT_bool, bg='#FAFBFC')
        self.isPT.grid(row=1, column=0, sticky='W', padx=10, pady=3)
        
        self.exp_name_text = Label(tab1, text='Name for Results Repertory', bg='#FAFBFC')
        self.exp_name_text.grid(row=2, column=0, sticky='W', padx=5, pady=10)
        self.exp_name = Entry(tab1)
        self.exp_name.grid(row=2, column=1, columnspan=10, sticky="WE", pady=3)
        self.exp_name.insert(0, "results_files_loc")
        
        # Select statistic to use for heatmaps
        self.stats_choice = Label(tab1, text = "Statistical value used for heatmaps", bg='#FAFBFC')
        self.stats_choice.grid(row=3, column=0, sticky="W", pady=3, ipadx=1, padx=5)
        self.list_stats_choice = ["median", "mean"]
        self.widget_stats_choice = ttk.Combobox(tab1, values=self.list_stats_choice, state="readonly")
        self.widget_stats_choice.grid(row=4, column=0, sticky="W", pady=3, ipadx=1, padx=5)
        self.widget_stats_choice.current(0)
        self.choice_stats = 'median'
        self.method_choice_stats = statistics.mean
        
        # Select if we want boxplot or histogram
        self.hist_box = Label(tab1, text = "Distribution Representation", bg='#FAFBFC')
        self.hist_box.grid(row=5, column=0, sticky="W", pady=3, ipadx=1, padx=5)
        self.list_hist_box = ["density histogram", "boxplot"]
        self.widget_hist_box = ttk.Combobox(tab1, values=self.list_hist_box, state="readonly")
        self.widget_hist_box.grid(row=6, column=0, sticky="W", pady=3, ipadx=1, padx=5)
        self.widget_hist_box.current(0)
        self.hist_box_choice = 'density histogram'
        self.use_boxplot = True

                
        # ------- LOCALIZATIONS ANALYSIS TAB -------
        self.run_exp_bool = BooleanVar()
        self.run_exp_bool.set(False)
        self.run_exp = Button(tab2, text='Run Density Measurement', command=self.do_run_loc_exp)
        self.run_exp.grid(row=0, column=0, sticky="WE", pady=3, ipadx=1, padx=5)
        
        self.run_exp_bool = BooleanVar()
        self.run_exp_bool.set(False)
        self.run_exp = Button(tab2, text='Get Cumulative Number of Clusters', command=self.do_run_cum_num_clus)
        self.run_exp.grid(row=1, column=0, sticky="WE", pady=3, ipadx=1, padx=5)

        self.run_exp_bool = BooleanVar()
        self.run_exp_bool.set(False)
        self.run_exp = Button(tab2, text='Get Photophysics Plots', command=self.do_photophysics_analysis)
        self.run_exp.grid(row=2, column=0, sticky="WE", pady=3, ipadx=1, padx=5)
        self.poca_files = None
        self.csv_intensity_files = None
        self.csv_frame_files = None
        self.csv_sigma_files = None      
        
          
        # ------- CLUSTERS ANALYSIS TAB -------
        self.get_experiment_heatmap_bool = BooleanVar()
        self.get_experiment_heatmap_bool.set(False)
        self.get_experiment_heatmap = Button(tab3, text='Get Experiment Heatmap', command=self.do_heatmap_whole_exp)
        self.get_experiment_heatmap.grid(row=0, column=0, sticky="WE", pady=3, ipadx=1, padx=5)

        self.get_one_heatmap_bool = BooleanVar()
        self.get_one_heatmap_bool.set(False)
        self.get_one_heatmap = Button(tab3, text='Get Photophysics Parameter Heatmap', command=self.do_one_heatmap)
        self.get_one_heatmap.grid(row=1, column=0, sticky="WE", pady=3, ipadx=1, padx=5)
        
        self.index_we_want = []
        self.phot_parameters = ['ON times', 'OFF times', "Intensity_loc", 'total ON',
                                'blinks', 'intensity', '# seq ON', '# seq OFF', 'Loc_Precision']

        self.options = ["Length ON times", "Length OFF times", "Intensity per Loc.", "Total ON time",
                "Num. Blinks", "Intensity per Clus.", "Num. ON time", "Num. OFF time", "Loc. Precision"]

        self.checkbox_vars = []
        self.checkboxs = list()
        for i in range(0,4):
            var = IntVar()
            self.checkbox_vars.append(var)
            checkbox = Checkbutton(tab3, text=self.options[i], variable=var, bg='#FAFBFC')
            checkbox.grid(row=i+2, column=0, sticky='W')
            self.checkboxs.append(checkbox)
        for i in range(4,8):
            var = IntVar()
            self.checkbox_vars.append(var)
            checkbox = Checkbutton(tab3, text=self.options[i], variable=var, bg='#FAFBFC')
            checkbox.grid(row=i-2, column=1, sticky='W')
            self.checkboxs.append(checkbox)
        # Loc precision
        var = IntVar()
        self.checkbox_vars.append(var)
        checkbox = Checkbutton(tab3, text=self.options[8], variable=var, bg='#FAFBFC')
        checkbox.grid(row=7, column=0, sticky='W')
        self.checkboxs.append(checkbox)
        
        self.check_everything = Button(tab3, text='Check Everything', command=self.select_all, bg='#FAFBFC')
        self.check_everything.grid(row=8, column=0, sticky='W')
        
        
        # ------- 96-WELL PLATE ANALYSIS TAB -------
        self.get_96heatmap_bool = BooleanVar()
        self.get_96heatmap_bool.set(False)
        self.get_96heatmap = Button(tab4, text='Get Photophysics Row Heatmaps', command=self.do_96heatmap)
        self.get_96heatmap.grid(row=1, column=0, sticky="WE", pady=3, ipadx=1, padx=5)

        self.get_one_96heatmap_bool = BooleanVar()
        self.get_one_96heatmap_bool.set(False)
        self.get_one_96heatmap = Button(tab4, text='Get Photophysics Parameter Heatmap', command=self.do_one_96heatmap)
        self.get_one_96heatmap.grid(row=1, column=1, sticky="WE", pady=3, ipadx=1, padx=5)
        
        self.index_we_want = []
        self.phot_parameters = ['ON times', 'OFF times', "Intensity_loc", 'total ON',
                                'blinks', 'intensity', '# seq ON', '# seq OFF', 'Loc_Precision']

        self.options = ["ON times Duration", "OFF times Duration", "Intensity per Loc.", "Total ON Time",
                "Number of Blinks", "Intensity per Molecule", "Number ON times", "Number OFF times", "Loc. Precision"]

        self.checkbox_vars96 = []
        self.checkboxs96 = list()
        for i in range(0,4):
            var = IntVar()
            self.checkbox_vars96.append(var)
            checkbox = Checkbutton(tab4, text=self.options[i], variable=var, bg='#FAFBFC')
            checkbox.grid(row=i+2, column=0, sticky='W')
            self.checkboxs96.append(checkbox)
        for i in range(4,8):
            var = IntVar()
            self.checkbox_vars96.append(var)
            checkbox = Checkbutton(tab4, text=self.options[i], variable=var, bg='#FAFBFC')
            checkbox.grid(row=i-2, column=1, sticky='W')
            self.checkboxs96.append(checkbox)
        # Loc precision
        var = IntVar()
        self.checkbox_vars96.append(var)
        checkbox = Checkbutton(tab4, text=self.options[8], variable=var, bg='#FAFBFC')
        checkbox.grid(row=7, column=0, sticky='W')
        self.checkboxs96.append(checkbox)
        self.check_everything96 = Button(tab4, text='Check Everything', command=self.select_all96, bg='#FAFBFC')
        self.check_everything96.grid(row=8, column=0, sticky='WE', pady=3, ipadx=1, padx=5)

        self.run_time_fov_bool = BooleanVar()
        self.run_time_fov_bool.set(False)
        self.run_time_fov = Button(tab4, text='Time Distribution for each well (FOVs)', command=self.plot_FOVs_well_time_distribution)
        self.run_time_fov.grid(row=9, column=0, sticky="WE", pady=3, ipadx=1, padx=5)

        self.run_time_row_bool = BooleanVar()
        self.run_time_row_bool.set(False)
        self.run_time_row = Button(tab4, text='Time Distribution for each row', command=self.plot_rows_time_distribution)
        self.run_time_row.grid(row=9, column=1, sticky="WE", pady=3, ipadx=1, padx=5)

        self.run_stat_test_fov_bool = BooleanVar()
        self.run_stat_test_fov_bool.set(False)
        self.run_stat_test_fov = Button(tab4, text='Run Stat Tests for FOVs', command=self.run_fov_stats_test)
        self.run_stat_test_fov.grid(row=10, column=0, sticky="WE", pady=3, ipadx=1, padx=5)

        self.run_stat_test_row_bool = BooleanVar()
        self.run_stat_test_row_bool.set(False)
        self.run_stat_test_row = Button(tab4, text='Run Stat Tests for rows', command=self.run_row_stats_test)
        self.run_stat_test_row.grid(row=10, column=1, sticky="WE", pady=3, ipadx=1, padx=5)   

        self.run_fps_heatmap_bool = BooleanVar()
        self.run_fps_heatmap_bool.set(False)
        self.run_fps_heatmap = Button(tab4, text='Get FPs Global Heatmap', command=self.do_fps_heatmap)
        self.run_fps_heatmap.grid(row=11, column=0, sticky="WE", pady=3, ipadx=1, padx=5)

    def load_molecule_data(self):
        """
        Load PALMTracer files (ending with locPALMTracer.txt), PoCA files (cleaned PT and csv files)
        """        
        try:
            self.repertory_path = askdirectory(title="Select your Source directory")
            self.index = [f for f in os.listdir(self.repertory_path+'/') if os.path.isdir(os.path.join(self.repertory_path+'/', f))]
            self.files = get_PALMTracer_files(self.repertory_path)
            self.poca_files = get_poca_files(self.repertory_path)
            self.csv_intensity_files = get_csv_poca_intensity_files(self.repertory_path)
            self.csv_frame_files = get_csv_poca_frame_files(self.repertory_path)
            self.csv_sigma_files = get_csv_poca_sigma_files(self.repertory_path)
            self.isPT_bool.set(True)
            print("Done:", "'"+self.repertory_path+"'", 'has been loaded')
        except:
            print("No directory selected")


    def do_run_loc_exp(self):
        """
        Run the localization data analysis and plot the density through time
        """        
        for i in self.files:
            if self.isPT_bool.get() == False:
                self.laser = self.index
            # Get laser and well name
            idx = os.sep.join(i.replace('/SR_001.MIA/locPALMTracer.txt', '').split('/')[-2:])
            # Get data about localisation, density...
            cum_loc_per_frame, density_per_frame, avg_density = do_analysis_for_one_acquisition(i)
            # Plot and saving plots
            save_loc_as_pdf(i, idx, density_per_frame, cum_loc_per_frame, avg_density, self.exp_name.get())
        print("Analysis Done!")
        
                
    def do_run_cum_num_clus(self):
        do_cumulative_number_clusters(self.poca_files, self.exp_name.get(), self.isPT_bool.get())
        print("Cumulative Clusters Analysis Done!")
    
    
    def do_photophysics_analysis(self):
        self.select_stats_method_visu()
        do_photophysics_parameters_plotting(self.poca_files, self.csv_frame_files, self.csv_intensity_files, self.csv_sigma_files, 
                                            self.exp_name.get(), self.isPT_bool.get(), boxplot=self.use_boxplot)
        print("Cluster Photophysics Plotting Done!")

        
    def do_heatmap_whole_exp(self):
        self.select_stats_method_heatmap()
        do_heatmap_photophysics_parameters(self.exp_name.get(), self.poca_files, self.csv_frame_files, self.csv_intensity_files, 
                                           self.csv_sigma_files, self.isPT_bool.get(), stats=self.method_choice_stats)
        print("Heatmap for the Whole Experiment Done!")

        
    def do_one_heatmap(self):
        self.select_stats_method_heatmap()
        self.checkbox_values = [option for option, var in zip(range(0, 9), self.checkbox_vars) if var.get()]
        self.index_we_want = [self.phot_parameters[i] for i in self.checkbox_values]
        do_heatmap_one_photophysics_parameter(self.exp_name.get(), self.index_we_want, self.poca_files, self.csv_frame_files, 
                                              self.csv_intensity_files, self.csv_sigma_files,self.isPT_bool.get(), 
                                              stats=self.method_choice_stats)
        print("Heatmap(s) for Selected Parameters Done!")


    def do_96heatmap(self):
        self.select_stats_method_heatmap()
        idx = list(map(str, range(1, 13)))  # Indices possibles, ajustables selon vos besoins
        cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))  # Colonnes possibles, ajustables
        # Identifier les colonnes disponibles dans les fichiers
        available_cols = set()
        for f in self.poca_files + self.csv_sigma_files + self.csv_intensity_files + self.csv_frame_files:
            for col in cols:
                if any(f"/{col}{i}/" in f for i in idx):
                    available_cols.add(col)
        # Convertir le set en liste ordonnée
        available_cols = sorted(list(available_cols))
        # Maintenant, bouclez uniquement sur les colonnes disponibles
        for col in available_cols:
            relevant_poca_files = [f for f in self.poca_files if any(f"/{col}{i}/" in f for i in idx)]
            relevant_sigma_files = [f for f in self.csv_sigma_files if any(f"/{col}{i}/" in f for i in idx)]
            relevant_intensity_files = [f for f in self.csv_intensity_files if any(f"/{col}{i}/" in f for i in idx)]
            relevant_frame_files = [f for f in self.csv_frame_files if any(f"/{col}{i}/" in f for i in idx)]
            # Appel de la fonction pour obtenir les données de la carte thermique
            do_96heatmap_photophysics_parameters(self.exp_name.get(), relevant_poca_files, relevant_frame_files, 
                                                 relevant_intensity_files, relevant_sigma_files, col, 
                                                 self.isPT_bool.get(), stats=self.method_choice_stats)
        print("HCS Heatmaps for all Rows Done!")



        
    def do_one_96heatmap(self):
        self.select_stats_method_heatmap()
        self.checkbox_values = [option for option, var in zip(range(0, 9), self.checkbox_vars96) if var.get()]
        self.index_we_want = [self.phot_parameters[i] for i in self.checkbox_values]
        do_96heatmap_one_photophysics_parameter(self.exp_name.get(), self.index_we_want, self.poca_files, self.csv_frame_files, 
                                              self.csv_intensity_files, self.csv_sigma_files,self.isPT_bool.get(), 
                                              stats=self.method_choice_stats)
        print("HCS Heatmap(s) for Selected Parameters Done!")


    def select_all(self):
        for i in self.checkboxs:
            i.invoke()
            
    def select_all96(self):
        for i in self.checkboxs:
            i.invoke()

    def select_all_super(self):
        for i in self.checkboxs_super:
            i.invoke()
                    
    def select_stats_method_heatmap(self):
        # Obtenir l'élément sélectionné
        self.choice_stats = self.widget_stats_choice.get()
        if self.choice_stats == 'mean':
            self.method_choice_stats = statistics.mean
        else:
            self.method_choice_stats = statistics.median


    def select_stats_method_visu(self):
        self.hist_box_choice = self.widget_hist_box.get()
        if self.hist_box_choice == 'boxplot':
            self.use_boxplot = True
        else:
            self.use_boxplot = False


    def plot_FOVs_well_time_distribution(self):
        idx = list(map(str, range(1, 13)))
        cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
        all_wells = fusion(cols, idx)
        for i in all_wells:
            if any(i in j for j in self.poca_files):
                time_distribution_same_well(self.repertory_path+'/'+i, i, self.exp_name.get())
                fov_photophysics_distribution(self.repertory_path+'/'+i, self.exp_name.get())
        print("Time Evolution between and for each well FOVs Done")


    def plot_rows_time_distribution(self):
        idx = list(map(str, range(1, 13)))  # Indices possibles, ajustables selon vos besoins
        cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))  # Colonnes possibles, ajustables
        for col in cols:
            super_table = []  # Contient les données pour toutes les `wells` du `col`
            wells = []
            relevant_files = [f for f in self.poca_files if any(f"/{col}{i}/" in f for i in idx)]
            # Grouper les fichiers pertinents par puit
            for well in [f"{col}{i}" for i in idx]:
                well_files = [f for f in relevant_files if f"/{well}/" in f]
                if well_files:
                    combined_data = pd.DataFrame()
                    for file in well_files:
                        raw_data = read_poca_files(file)
                        filtered_data = raw_data[raw_data['total ON'] < np.quantile(raw_data['total ON'], 1)]
                        filtered_data = filtered_data[filtered_data['blinks'] < np.quantile(filtered_data['blinks'], 1)]
                        filtered_data["avg_on"] = np.array(filtered_data['total ON'] / filtered_data['# seq ON'])
                        filtered_data["photon_loc"] = np.array(filtered_data['intensity'] / filtered_data['total ON'])
                        grouped_data = filtered_data.groupby(pd.cut(filtered_data['frame'], list(range(0, 5001, 500)))).mean()
                        combined_data = pd.concat([combined_data, grouped_data], ignore_index=True)
                    # well_photophysics_distribution(combined_data, well, self.exp_name.get())
                    super_table.append(combined_data)
                    wells.append(well)
            if super_table:
                time_distribution_same_row(super_table, self.exp_name.get(), col, wells)
        print("Time Evolution between and for Plate Rows Done!")


    def run_fov_stats_test(self):
        idx = list(map(str, range(1, 13)))
        cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
        all_wells = fusion(cols, idx)
        self.checkbox_values = [option for option, var in zip(range(0, 9), self.checkbox_vars96) if var.get()]
        self.index_we_want = [self.phot_parameters[i] for i in self.checkbox_values]
        for i in all_wells:
            if any(i in j for j in self.poca_files):
                poca_files = get_poca_files(self.repertory_path+'/'+i)
                csv_intensity_files = get_csv_poca_intensity_files(self.repertory_path+'/'+i)
                csv_frame_files = get_csv_poca_frame_files(self.repertory_path+'/'+i)
                csv_sigma_files = get_csv_poca_sigma_files(self.repertory_path+'/'+i)
                do_stat_tests_fovs(self.exp_name.get(), i, self.index_we_want, poca_files, csv_frame_files, 
                             csv_intensity_files, csv_sigma_files)
        print("Statistical Tests for FOVs of each well done!")


    def run_row_stats_test(self):
        idx = list(map(str, range(1, 13)))
        cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
        self.checkbox_values = [option for option, var in zip(range(0, 9), self.checkbox_vars96) if var.get()]
        self.index_we_want = [self.phot_parameters[i] for i in self.checkbox_values]
        for col in cols:
            wells_data = {}
            for i in idx:
                well = f"{col}{i}"
                if any(f"/{well}/" in j for j in self.poca_files):                    
                    # Récupérer les fichiers pour le puit
                    poca_files = get_poca_files(self.repertory_path + '/' + well)
                    csv_intensity_files = get_csv_poca_intensity_files(self.repertory_path + '/' + well)
                    csv_frame_files = get_csv_poca_frame_files(self.repertory_path + '/' + well)
                    csv_sigma_files = get_csv_poca_sigma_files(self.repertory_path + '/' + well)
                    # Concaténer les données de chaque type
                    well_poca_data = pd.DataFrame()
                    for file in poca_files:
                        raw_data = read_poca_files(file)
                        well_poca_data = pd.concat([well_poca_data, raw_data], ignore_index=True)
                    wells_data[well] = {
                        "poca_data": well_poca_data,
                        "frame_data": csv_frame_files,
                        "intensity_data": csv_intensity_files,
                        "sigma_data": csv_sigma_files
                    }
            if wells_data:
                # Appliquer les tests statistiques sur les données regroupées par puit
                do_stats_tests_row(self.exp_name.get(), col, self.index_we_want, wells_data)
        print("Statistical Tests for FOVs of each well and between wells in the same row done!")



    def do_fps_heatmap(self):
        self.select_stats_method_heatmap()
        # Liste des protéines
        protein_list = ["mEos4b-L93M", "mEos3.2", "mEosEM", "pcSTAR", "Dendra2", "mAple3", "mEos4b"]
        all_protein_data = {}

        for protein in protein_list:
            relevant_poca_files = [f for f in self.poca_files if f"/{protein}/" in f]
            relevant_sigma_files = [f for f in self.csv_sigma_files if f"/{protein}/" in f]
            relevant_intensity_files = [f for f in self.csv_intensity_files if f"/{protein}/" in f]
            relevant_frame_files = [f for f in self.csv_frame_files if f"/{protein}/" in f]

            if relevant_poca_files:
                protein_data = do_fps_heatmap(
                    relevant_poca_files, relevant_frame_files,
                    relevant_intensity_files, relevant_sigma_files,
                    self.isPT_bool.get(), stats=self.method_choice_stats
                )
                all_protein_data[protein] = protein_data

        generate_fps_heatmap(all_protein_data, self.exp_name.get())
        print("FPs Heatmap for all proteins done!")