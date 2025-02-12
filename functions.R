# Define the write_fasta function
write_fasta <- function(target_allele, df, id_col, seq_col, file_path) {
  # Open a file connection
  con <- file(file_path, "w")
  
  # Create FASTA entries
  fasta_entries <- df |>
    filter(str_detect(MHC_Allele, target_allele)) |>  
    group_by(Epitope) |>
    slice_head(n = 1) |>
    ungroup() |>
    rowwise() |>  # Ensure row-wise operation
    mutate(fasta_entry = paste0(">", !!sym(id_col), "\n", !!sym(seq_col))) |> 
    pull(fasta_entry)
  
  # Write to file
  writeLines(fasta_entries, con)
  
  # Close the file connection
  close(con)
}

cross_reactive_test <- function(percent_identity, autoimmune_data, pathogen_df, target_allele){
  # Running similarity check with % acceptable errors on our now filtered proteome data
  cross_reactive_results <- autoimmune_data |>
    mutate(matches = map(Epitope, ~{
      peptide_length <- nchar(.x)  # Get peptide length
      max_mismatch <- round(peptide_length * (1 - percent_identity/100))  # Compute max mismatches
      
      matched_indices <- which(vcountPattern(pattern = .x, 
                                             subject = pathogen_df$sequence,
                                             max.mismatch = max_mismatch,
                                             with.indels = TRUE) > 0)
      list(Pathogen_ID = pathogen_df$header[matched_indices],
           count = length(matched_indices))
    })) |>
    unnest_wider(matches) |>
    filter(count >= 1) |>
    arrange(desc(count)) |>
    mutate(Pathogen_ID = map(Pathogen_ID, unique)) |>
    unnest(Pathogen_ID) |>
    mutate(Pathogen_ID = as.character(Pathogen_ID))
  
  return(cross_reactive_results)
}

MHC_motif_finder <- function(MHC_allele_pattern){
  
  mhc_ref_list <- read_delim("./Data/all_logos/pseudo_mhc_list", delim = " ", col_names = c("filename", "MHC_allele"))
  
  logo_filename <- mhc_ref_list |>
    filter(str_detect(MHC_allele, MHC_allele_pattern)) |>
    pull("filename")
  
  pdf_location <- paste("./Data/all_logos/", logo_filename, "_el.pdf", sep = "")
  
  # Convert first page of PDF to an image
  pdf_img <- image_read_pdf(pdf_location)
  
  return(pdf_img)
}
