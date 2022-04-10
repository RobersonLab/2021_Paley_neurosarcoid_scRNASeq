# title: "paley_neurosarcoid_2020_08 - Prefilter to find Doublets"
# author: "Elisha Roberson"
# date: "2020-08-27"

# Analysis related to - paley_neurosarcoid_2020_08

# packages
library( here )
library( tidyverse )
library( Seurat )
library( DoubletFinder )

# source
source( file = here::here( "src", "shared_r_functions.R" ) )

# decide graphics parameter based on OS
graphics_device_type <- case_when(
  .Platform$OS.type == "windows" ~ "windows",
  TRUE ~ "cairo"
)

image_width <- 1280
image_height <- 720

# read project info
project_info <- read_tsv( file = here( "info", "unified_project_info.txt" ) )

# double finder loop
for ( row_idx in 1:nrow( project_info ) ) {
  rgid <- project_info %>%
    pull( rgid ) %>%
    .[ row_idx ]

  output_filename <- paste0( rgid, "_doublets.tsv.gz" )

  # load data
  count_data <- Read10X( data.dir = here( "data", "raw_umi", rgid ), strip.suffix = TRUE )

  count_meta <- data.frame( cell_name = colnames( count_data ), Status = "mock" ) %>%
    column_to_rownames( "cell_name" )

  # convert to seurat
  seurat_data <- CreateSeuratObject( counts = count_data,
                                     project = "Null",
                                     min.cells = min_cell_cutoff,
                                     min.features = min_feature_cutoff,
                                     meta.data = count_meta )

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

  seurat_data <- subset( seurat_data, nFeature_RNA > min_feature_cutoff & percent.mt < mt_cutoff_value )

  # basic transforms
  seurat_data <- SCTransform( seurat_data )
  seurat_data <- RunPCA( seurat_data )
  seurat_data <- RunUMAP( seurat_data, dims = 1:10 )

  # parameter sweep
  sweep_res_list <- paramSweep_v3( seurat_data, PCs = 1:10, sct = TRUE )

  # pK identification
  sweep_stats <- summarizeSweep( sweep_res_list, GT = FALSE )
  found_pK <- find.pK( sweep_stats )

  pK_value <- arrange( found_pK, desc( BCmetric ) ) %>%
    pull( pK ) %>%
    .[ 1 ] %>%
    as.character( . ) %>%
    as.numeric( . )

  # number of expected doublets
  nExp_poi <- round( expected_doublet_fraction * ncol( seurat_data ) )

  # doublet finder with calculated pK and expected doublets
  seurat_data <- doubletFinder_v3( seurat_data,
                                   PCs = 1:10,
                                   pN = 0.25,
                                   pK = pK_value,
                                   nExp = nExp_poi,
                                   reuse.pANN = FALSE,
                                   sct = TRUE )

  # column names of DoubletFinder vary by the parameters
  # this step just finds them for us
  pANN_col_num <- seurat_data@meta.data %>%
    colnames( . ) %>%
    str_detect( string = ., pattern = "^pANN" ) %>%
    which( . )

  pANN_col_name <- seurat_data@meta.data %>%
    colnames( . ) %>%
    .[ pANN_col_num ]

  DF_col_num <- seurat_data@meta.data %>%
    colnames( . ) %>%
    str_detect( string = ., pattern = "^DF.classifications" ) %>%
    which( . )

  DF_col_name <- seurat_data@meta.data %>%
    colnames( . ) %>%
    .[ DF_col_num ]

  # show the plots of doublets
  doublet_plot <- DimPlot( seurat_data,
                           reduction = "umap",
                           group.by = DF_col_name,
                           pt.size = 2 ) + ggtitle( rgid )

  jpeg( width = image_width,
        height = image_height,
        type = graphics_device_type,
        filename = here( "results", "figures", "doublets" , paste0( rgid, "_doublet_plot.jpeg" ) ) )

  print( doublet_plot )

  dev.off()

  # pull some of this data into another dataframe
  my_df_doublets <- data.frame( classification = seurat_data@meta.data[ , DF_col_name ],
                                pANN = seurat_data@meta.data[ , pANN_col_num ],
                                cell_id = rownames( seurat_data@meta.data ) )
  rownames( my_df_doublets ) = my_df_doublets$cell_id

  # get names of potential doublets
  name_of_bad_cells <- my_df_doublets %>%
    filter( classification == "Doublet" ) %>%
    pull( cell_id ) %>%
    as.character( . ) %>%
    data.frame( cell_id = ., category = "Doublet" ) %>%
    write_tsv( x = ., path = here( "results", "doublets", output_filename ) )

  # heatmap the doublets to look for overlap
  heatmap_genelist <- c( 'CD3D', 'CD3E', 'CD19', 'CD79A', 'LYZ', 'NCR1', 'CD68' )

  doublet_list <- my_df_doublets %>%
    filter( classification == "Doublet" ) %>%
    pull( cell_id )

  heatmap_plot <- DoHeatmap( seurat_data,
                             features = heatmap_genelist,
                             cells = doublet_list )

  jpeg( width = image_width,
        height = image_height,
        type = graphics_device_type,
        filename = here( "results", "figures", "doublets" , paste0( rgid, "_doublet_heatmap.jpeg" ) ) )

  print( heatmap_plot )

  dev.off()
}

# session info
Sys.time()

getwd()

sessionInfo()
