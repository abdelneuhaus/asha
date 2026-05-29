#' Read a locPALMTracer file
#'
#' Reads a tab-separated locPALMTracer file, skips the first two header rows, 
#' and enforces integer types for specific columns.
#'
#' @param file Character string. Path to the locPALMTracer text file.
#' @return A data.frame with correctly typed columns.
#' @export
read_locPALMTracer_file <- function(file) {
  data <- read.delim(file, sep = '\t', skip = 2, stringsAsFactors = FALSE)
  
  data$id <- as.integer(data$id)
  data$Plane <- as.integer(data$Plane)
  data$Index <- as.integer(data$Index)
  
  return(data)
}

#' Read and clean a PoCA file
#'
#' Reads a PoCA CSV file, replaces outliers in the `# seq OFF` column with 
#' values from the `blinks` column, and removes the last filler column.
#'
#' @param file Character string. Path to the PoCA CSV file.
#' @return A cleaned data.frame.
#' @export
read_poca_files <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  if ("# seq OFF" %in% colnames(df) && "blinks" %in% colnames(df)) {
    condition <- df[["# seq OFF"]] > 5000
    condition[is.na(condition)] <- FALSE
    df[["# seq OFF"]][condition] <- df[["blinks"]][condition]
  }
  
  df <- df[, -ncol(df), drop = FALSE]
  
  return(df)
}

#' Read a generic PoCA CSV into a list of vectors
#'
#' @param file Character string. Path to the CSV file.
#' @return A list where each element represents a row from the data.frame.
#' @export
read_csv_poca <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  return(split(as.matrix(df), row(df)))
}

#' Extract the frame column from a PoCA file
#'
#' @param file Character string. Path to the PoCA CSV file.
#' @return A list containing the values of the 'frame' column.
#' @export
read_cluster_from_PoCA <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  df <- df[, -ncol(df), drop = FALSE]
  
  return(as.list(df[["frame"]]))
}


#' Recursively find PoCA merged or cleaned files
#'
#' Scans a directory for specific PoCA output files. Filters results to only 
#' include files located inside folders containing "MIA" in their name.
#'
#' @param repertory Character string. The root directory to search.
#' @param new Logical. If TRUE, searches for 'locPALMTracer_merged.txt'. 
#'   If FALSE, searches for 'cleaned.txt'. Default is TRUE.
#' @return A character vector of normalized file paths.
#' @export
get_poca_files <- function(repertory, new = TRUE) {
  suffix <- if (new) "locPALMTracer_merged\\.txt$" else "cleaned\\.txt$"
  all_files <- list.files(repertory, pattern = suffix, recursive = TRUE, full.names = TRUE)
  
  valid_files <- all_files[grepl("MIA", dirname(all_files))]
  valid_files <- gsub("\\\\", "/", valid_files)
  
  return(valid_files)
}


#' Find intensity CSV files
#'
#' @param repertory Character string. Root directory to search.
#' @return A character vector of file paths.
#' @export
get_csv_poca_intensity_files <- function(repertory) {
  files <- list.files(repertory, pattern = "intensity\\.csv$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}


#' Find sigmaXY CSV files
#'
#' @param repertory Character string. Root directory to search.
#' @return A character vector of file paths.
#' @export
get_csv_poca_sigma_files <- function(repertory) {
  files <- list.files(repertory, pattern = "sigmaXY\\.csv$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}


#' Find frame CSV files
#'
#' @param repertory Character string. Root directory to search.
#' @return A character vector of file paths.
#' @export
get_csv_poca_frame_files <- function(repertory) {
  files <- list.files(repertory, pattern = "frame\\.csv$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}

#' Find locPALMTracer files
#'
#' @param repertory Character string. Root directory to search.
#' @return A character vector of file paths.
#' @export
get_PALMTracer_files <- function(repertory) {
  files <- list.files(repertory, pattern = "locPALMTracer\\.txt$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}


#' Find trcPALMTracer files
#'
#' @param repertory Character string. Root directory to search.
#' @return A character vector of file paths.
#' @export
get_trcPALMTracer_files <- function(repertory) {
  files <- list.files(repertory, pattern = "trcPALMTracer\\.txt$", recursive = TRUE, full.names = TRUE)
  return(gsub("\\\\", "/", files))
}