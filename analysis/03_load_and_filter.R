# title: "paley_neurosarcoid_2020_08 - Load and filter cells"
# author: "Elisha Roberson"
# date: "2020-08-27"

# Analysis related to - paley_neurosarcoid_2020_08

# load libraries
library( here )
library( tidyverse )
library( Seurat )

# Source
source( file = here( "src", "shared_r_functions.R" ) )

# project info
project_info <- read_tsv( file = here::here( "info", "unified_project_info.txt" ) )

# Read raw data
# Filter features / mito
# Remove doublets

seurat_object_name_list = character( 0 )
seurat_cell_id_list = character( 0 )

cell_numbers = data.frame( rgsm = project_info$rgsm, raw_cells = 0, filtered_cells = 0 )
rownames( cell_numbers ) = cell_numbers$rgsm

for ( row_idx in 1:nrow( project_info ) ) {
  # current sample info
  rgid <- project_info %>%
    pull( rgid ) %>%
    .[ row_idx ]

  rgsm <- project_info %>%
    pull( rgsm ) %>%
    .[ row_idx ]

  gen_status <- project_info %>%
    pull( gen_status ) %>%
    .[ row_idx ]

  status <- project_info %>%
    pull( status ) %>%
    .[ row_idx ]

  tissue <- project_info %>%
    pull( tissue ) %>%
    .[ row_idx ]

  diagnosis <- project_info %>%
    pull( diagnosis ) %>%
    .[ row_idx ]

  doublet_filename <- paste0( rgid, "_doublets.tsv.gz" )

  # load doublet list
  doublet_cells <- read_tsv( file = here( "results", "doublets", doublet_filename ) ) %>%
    pull( cell_id )

  # load data
  count_data <- Read10X( data.dir = here( "data", "raw_umi", rgid ), strip.suffix = TRUE )

  tcr_meta <- read_tsv( file = here( 'results', 'tcr_meta_data', paste0( rgid, '_tcr_qc.tsv.gz' ) ) ) %>%
    dplyr::rename( cell_name = barcode )

  count_meta <- data.frame( cell_name = colnames( count_data ), Status = status, Tissue = tissue, GeneralStatus = gen_status, Diagnosis = diagnosis, Sample = rgsm ) %>%
  merge( x = ., y = tcr_meta, all.x = TRUE, by = "cell_name" ) %>%
  column_to_rownames( "cell_name" )

  # convert to seurat
  seurat_data <- CreateSeuratObject( counts = count_data,
                                     project = project_name,
                                     min.cells = min_cell_cutoff,
                                     min.features = min_feature_cutoff,
                                     meta.data = count_meta )
  # raw counts
  cell_numbers[ rgsm, "raw_cells" ] = length( colnames( seurat_data ) )

  # remove intermediate objects
  rm( count_data )
  rm( count_meta )

  # Basic feature / mitochondrial QC
  gene_name_list <- rownames( seurat_data )

  mt_string <- case_when(
    length( which( str_detect( string = gene_name_list, pattern = "^Mt-" ) ) ) > 0 ~ "^Mt-",
    length( which( str_detect( string = gene_name_list, pattern = "^MT-" ) ) ) > 0 ~ "^MT-",
    TRUE ~ NULL
  )

  if ( is.null( mt_string ) ) {
    stop( paste0( "No mitochondrial string detected" ) )
  }

  seurat_data[[ "percent.mt" ]] <- PercentageFeatureSet( seurat_data, pattern = mt_string )

  mt_cutoff_value <- case_when(
    length( which( seurat_data[[ "percent.mt" ]][,1] > 1.0 ) ) > 0 ~ 20.0,
    TRUE ~ 0.20
  )

  seurat_data <- subset( seurat_data, snFeature_RNA > min_feature_cutoff & percent.mt < mt_cutoff_value )

  # Remove doublets
  seurat_data <- seurat_data[ , setdiff( colnames( seurat_data ), doublet_cells ) ]

  # count filtered
  cell_numbers[ rgsm, "filtered_cells" ] = length( colnames( seurat_data ) )

  # save objects and move on
  object_name <- paste0( 'seurat_', row_idx )

  seurat_object_name_list[ row_idx ] = object_name
  seurat_cell_id_list[ row_idx ] = rgsm

  assign( x = object_name, value = seurat_data )
  rm( seurat_data )
}

# merge into seurat
seurat_object_list <- list()
for ( idx in 2:length( seurat_object_name_list ) ) {
  seurat_object_list[ idx - 1 ] <- base::get( seurat_object_name_list[ idx ] )
}

merged_raw <- merge( x = base::get( seurat_object_name_list[1] ),
                     y = seurat_object_list,
                     add.cell.ids = seurat_cell_id_list,
                     project = project_name )

# save R data objects
save( merged_raw, file = here( "results", "objects", "merged_raw.RData" ) )

# write cell numbers
write_tsv( x = cell_numbers, path = here( "results", "cell_counts.tsv" ) )

# Session info
Sys.time()

getwd()

sessionInfo()
