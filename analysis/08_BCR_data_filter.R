# title: bcr data filter
# author: Michael Paley

# load libaries
library( here )
library( ggplot2 )
library( cowplot )
library( dplyr )
library( tidyr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
library( readr )

#############
# constants #
#############
image_width <- 1280
image_height <- 720

graphics_device_type <- case_when(
  .Platform$OS.type == "windows" ~ "windows",
  TRUE ~ "cairo"
)

# Set core variables and read in "metadata" file
MetaData <- read_tsv( file = here( 'info', 'unified_project_info.tsv' ) )
nfiles <- nrow( MetaData )

# Load BCR data
BCRfiles <- dir( path = here( 'data', 'bcr' ) )
nBCRfiles <- length(BCRfiles)
BCR.samples <- BCRfiles

BCR_contig_list <- list()
for (i in 1:nBCRfiles) {
    BCR_contig_list[[i]] <- read_csv( file = here( 'data', 'bcr', BCRfiles[i], 'filtered_contig_annotations.csv' ) )
}

##################################################
# filter BCR VDJ contigs                         #
##################################################

input_BCR <- list()
for (i in 1:nBCRfiles) {
  input_BCR[[i]] <- BCR_contig_list[[i]] %>%
    dplyr::filter( full_length == TRUE ) %>%
    dplyr::filter( productive == TRUE ) %>%
    dplyr::filter( chain %in% c( 'IGH', 'IGK', 'IGL' ) ) %>%
    mutate( barcode = str_replace( string = barcode, pattern = "\\-1$", replacement = "" ) ) %>%
    mutate( contig_id = str_replace( string = contig_id, pattern = "\\-1", replacement =  "" ) )
}

bcr_summary <- list()
for (i in 1:nBCRfiles) {
  bcr_summary[[i]] <- dplyr::group_by( input_BCR[[i]], barcode, chain ) %>%
    dplyr::summarise( count = length( chain ) ) %>%
    pivot_wider( id_cols = barcode, values_from = count, names_from = chain, values_fill = 0 ) %>%
    select( barcode, IGH, IGL, IGK )
}

bcr_table <- list()
total <- list()
for (i in 1:nBCRfiles) {
  bcr_table[[i]] <- with( bcr_summary[[i]], table( IGH, IGL, IGK ) ) %>%
    as.data.frame( . )
  total[[i]] <- sum( bcr_table[[i]]$Freq )
  bcr_table[[i]] <- mutate( bcr_table[[i]], Percent = round( ( bcr_table[[i]]$Freq / total[[i]] ) * 100.0, 1 ) )
}

######### Plot QC of BCRs ##########

for (i in 1:nBCRfiles) {
  ggplot_bcr_freq <- ggplot( data = bcr_table[[i]], mapping = aes( IGH, IGL, IGK ) ) +
    theme_bw() +
    geom_tile( aes( fill = Freq ), colour = 'white' ) +
    scale_fill_gradient( low = "white", high = "red" ) +
    geom_text( aes( label = Freq ), size = 8 )

  jpeg( width = image_width,
        height = image_height,
        type = graphics_device_type,
        filename = here( 'results', 'figures', 'vdj_processing', paste0( BCR.samples[i], "_bcr_freq.jpeg" ) ) )

  print( ggplot_bcr_freq )

  dev.off()
}

rm(ggplot_bcr_freq)

for (i in 1:nBCRfiles) {
  ggplot_bcr_percent <- ggplot( data = bcr_table[[i]], mapping = aes( IGH, IGL, IGK ) ) +
    theme_bw() +
    geom_tile( aes( fill = Percent ), colour = 'white' ) +
    scale_fill_gradient( low = "white", high = "red" ) +
    geom_text( aes( label = Percent ), size = 8 )

  jpeg( width = image_width,
        height = image_height,
        type = graphics_device_type,
        filename = here( 'results', 'figures', 'vdj_processing', paste0( BCR.samples[i], "_bcr_percent.jpeg" ) ) )

  print( ggplot_bcr_percent )

  dev.off()
}

rm(ggplot_bcr_percent)

######### Plot QC of BCRs ##########

# Add summary of counts to big object
bcr_big_summary <- list()
for (i in 1:nBCRfiles) {
  bcr_table[[i]] <- mutate( bcr_table[[i]], SAMPLE = BCR.samples[i] ) # switched RGID/rgid to SAMPLE/Sample
  bcr_big_summary <- rbind( bcr_big_summary, bcr_table[[i]] )
}

# extract barcodes of good vdj & bad vdj
bcr_good <- list() # This will be a vector of barcodes for BCRs that pass the filter
bcr_bad <- list() # This will be a vector of barcodes for BCRs that DO NOT pass the filter
for (i in 1:nBCRfiles) {
  bcr_good[[i]] <- bcr_summary[[i]] %>%
    filter( IGH == 1 ) %>%
    filter( IGK + IGL == 1 ) %>%
    pull( barcode )

  bcr_bad[[i]] <- input_BCR[[i]] %>%
    pull( barcode ) %>%
    setdiff( x = ., y = bcr_good[[i]] )
}

number_bcrs <- list() # Creates a data frame of contig_id and the "chain_number"
for (i in 1:nBCRfiles) {
  number_bcrs[[i]] <- select( input_BCR[[i]], barcode, chain, contig_id ) %>%
    dplyr::group_by( barcode, chain ) %>%
    dplyr::summarise( chain_number = 1:length( chain ),
                      contig_id = contig_id ) %>%
    as.data.frame( . ) %>%
    select( contig_id, chain_number )
}

processed_BCR <- list() # Creates a data frame of with cell barcodes to add to Seurat object
for (i in 1:nBCRfiles) {
processed_BCR[[i]] <- select( input_BCR[[i]], contig_id, barcode, chain, raw_clonotype_id, raw_consensus_id, v_gene, d_gene, j_gene, cdr3, cdr3_nt ) %>%
  filter( barcode %in% bcr_good[[i]] ) %>%
  merge( x = number_bcrs[[i]], y = ., by = 'contig_id' ) %>%
  pivot_longer( cols = c( 'raw_consensus_id', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt' ) ) %>%
  mutate( name = paste0( chain, chain_number, "_", name ) ) %>%
  select( -contig_id, -chain, -chain_number ) %>%
  pivot_wider( names_from = name, values_from = value ) %>%
  merge( x = bcr_summary[[i]], y = ., by = 'barcode' ) %>%
  dplyr::rename( B_clonotype_id = raw_clonotype_id )
}

bcr_qc <- list()
for (i in 1:nBCRfiles) { # creates a .tsv files with each barcode (with '-1' at the tail) and whether it passed or failed the filter
  bcr_qc[[i]] <- data.frame( barcode = bcr_good[[i]], BCR_QC = "Pass" ) %>%
    rbind( ., data.frame( barcode = bcr_bad[[i]], BCR_QC = "Fail" ) )
  write_tsv( x = bcr_qc[[i]],
             here( 'results', 'bcr_qc', paste0( BCR.samples[i], '_bcr_qc.tsv.gz' ) ) )
}

for (i in 1:nBCRfiles) { # write processed_BCR file(s)
  processed_BCR[[i]] %>%
  write_tsv( x = .,
             file = here( 'results', 'vdj_meta_data', paste0( BCR.samples[i], '_bcr_qc.tsv.gz' ) ) )
}

BCR_frame <- data.frame()
for (i in 1:nBCRfiles) { # write processed_BCR file(s)
  BCR_frame <- mutate( processed_BCR[[i]], SAMPLE = BCR.samples[i] ) %>%
    dplyr::relocate( SAMPLE, .before = barcode ) %>%
    rbind( BCR_frame, . )
}
# write count summaries
write_tsv( x = bcr_big_summary, file = here( 'results', 'bcr_qc', 'global_bcr_count.tsv.gz' ) )

# write global vdj frame
write_tsv( x = BCR_frame, file = here( 'results', 'vdj_meta_data', 'global_bcr_processed.tsv.gz' ) )

# Remove intermediate objects
rm( BCR_contig_list )

# Table BCR clonotype
bcrfreq <- table( BCR_frame$B_clonotype_id , BCR_frame$SAMPLE ) %>%
  as.data.frame() %>%
  pivot_wider( . , names_from = Var2, values_from = Freq) %>%
  as.data.frame()

# Save calulations
write.csv( bcrfreq,
           file = here( 'results', 'calculations', 'VDJ', 'B.clonotype.counts.csv' ),
           quote = FALSE,
           row.names = FALSE )

######################### END ###############################

# session info
Sys.time()
getwd()
sessionInfo()

