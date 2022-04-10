#
#title: "paley_neurosarcoid_2020_08 - Integrate datasets"
#author: "Elisha Roberson"
#date: "2020-08-27"


# Analysis related to - paley_neurosarcoid_2020_08

# libraries
library( here )
library( tidyverse )
library( Seurat )

# source
source( file = here( "src", "shared_r_functions.R" ) )

# load raw
load( file = here( "results", "objects", "merged_raw.RData" ) )

# load info & cell counts
# get list of reference samples
info <- read_tsv( file = here( "info", "unified_project_info.txt" ) )
cell_counts <- read_tsv( file = here( "results", "cell_counts.tsv" ) )

get_top_rgsm <- function( df ) {
  top_n( x = df, n = 1, wt = filtered_cells ) %>%
    pull( rgsm ) %>%
    return( . )
}

top_samples <- merge( x = cell_counts, y = select( info, rgsm, tissue ), by = "rgsm" ) %>%
  group_by( tissue ) %>%
  nest( . ) %>%
  mutate( value = map_chr( data, get_top_rgsm ) ) %>%
  select( tissue, value )

reference_sample_names <- pull( top_samples, value )

## Find anchors
sc_list <- SplitObject( object = merged_raw, split.by = "Sample" )

for ( idx in 1:length( sc_list ) ) {
  sc_list[[ idx ]] <- SCTransform( sc_list[[ idx ]], vars.to.regress = "percent.mt", verbose = FALSE )
}

data_features <- SelectIntegrationFeatures( object.list = sc_list, nfeatures = 3000 )
sc_list <- PrepSCTIntegration( object.list = sc_list, anchor.features = data_features )

dataset_anchors <- FindIntegrationAnchors( object.list = sc_list, normalization.method = "SCT", anchor.features = data_features, reference = which( names( sc_list ) %in% reference_sample_names ) )

# integrate
integrated_data <- IntegrateData( anchorset = dataset_anchors, normalization.method = "SCT" )

save( integrated_data, file = here( "results", "objects", "integrated_data.RData" ) )

# Session info
Sys.time()

getwd()

sessionInfo()
