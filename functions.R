############################################################################################################################
# Function that takes in a MHC allele and outputs the logo plot for it

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

############################################################################################################################
# Filters sequences based on a given MHC allele anchor points, 
# and finds the binding core along with possible site of TCR binding

Allele_filterer <- function(df, allele_pattern, TCR_bind_pos){
  
  df_filtered <- df |>
    filter(str_detect(Seq, allele_pattern)) |>
    rowwise() |>
    mutate(start_index = str_locate(Seq, allele_pattern)[1]) |>
    mutate(IEDB_core = substring(Seq,start_index,start_index+8)) |>
    mutate(IEDB_tcore = paste0(strsplit(IEDB_core, "")[[1]][TCR_bind_pos],
                                  collapse = "")) |>
    ungroup()
}