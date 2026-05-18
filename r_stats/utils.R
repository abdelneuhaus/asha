# --- Fonctions de Lecture ---

read_locPALMTracer_file <- function(file) {
  # Lecture du fichier avec tabulation, en sautant les 2 premières lignes
  data <- read.delim(file, sep = '\t', skip = 2, stringsAsFactors = FALSE)
  
  # Conversion en entiers
  data$id <- as.integer(data$id)
  data$Plane <- as.integer(data$Plane)
  data$Index <- as.integer(data$Index)
  
  return(data)
}

read_poca_files <- function(file) {
  # check.names = FALSE permet de garder les noms de colonnes exacts comme "# seq OFF" 
  # au lieu de les transformer en "X..seq.OFF"
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Remplacement conditionnel : si '# seq OFF' > 5000, on prend la valeur de 'blinks'
  if ("# seq OFF" %in% colnames(df) && "blinks" %in% colnames(df)) {
    condition <- df[["# seq OFF"]] > 5000
    condition[is.na(condition)] <- FALSE # Sécurité contre les NA
    df[["# seq OFF"]][condition] <- df[["blinks"]][condition]
  }
  
  # df.iloc[:,:-1] -> On retourne toutes les colonnes sauf la dernière
  df <- df[, -ncol(df), drop = FALSE]
  
  return(df)
}

read_csv_poca <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  # L'équivalent R d'une liste de listes (values.tolist()) n'est pas très idiomatique. 
  # On le convertit en liste de vecteurs (une liste = une ligne) au cas où ton code 
  # en ait strictement besoin. Sinon, renvoyer simplement `df` est souvent mieux en R.
  return(split(as.matrix(df), row(df)))
}

read_cluster_from_PoCA <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Retire la dernière colonne
  df <- df[, -ncol(df), drop = FALSE]
  
  # Retourne la colonne 'frame' sous forme de liste (ou vecteur)
  return(as.list(df[["frame"]]))
}

# --- Fonctions de Recherche de Fichiers (Equivalent de os.walk) ---

get_poca_files <- function(repertory, new = TRUE) {
  # Définition du suffixe
  suffix <- if (new) "locPALMTracer_merged\\.txt$" else "cleaned\\.txt$"
  
  # list.files avec recursive = TRUE remplace os.walk
  all_files <- list.files(repertory, pattern = suffix, recursive = TRUE, full.names = TRUE)
  
  # On garde uniquement ceux qui ont "MIA" dans le chemin du dossier
  valid_files <- all_files[grepl("MIA", dirname(all_files))]
  
  # Normalisation des slash (R gère nativement les '/', on remplace les '\' éventuels)
  valid_files <- gsub("\\\\", "/", valid_files)
  
  return(valid_files)
}

get_csv_poca_intensity_files <- function(repertory) {
  files <- list.files(repertory, pattern = "intensity\\.csv$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}

get_csv_poca_sigma_files <- function(repertory) {
  files <- list.files(repertory, pattern = "sigmaXY\\.csv$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}

get_csv_poca_frame_files <- function(repertory) {
  files <- list.files(repertory, pattern = "frame\\.csv$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}

get_PALMTracer_files <- function(repertory) {
  files <- list.files(repertory, pattern = "locPALMTracer\\.txt$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}

get_trcPALMTracer_files <- function(repertory) {
  files <- list.files(repertory, pattern = "trcPALMTracer\\.txt$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}