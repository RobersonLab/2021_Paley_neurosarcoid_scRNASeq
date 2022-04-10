# title: "paley_neurosarcoid_2020_08 - Parse TCR info"
# author: "Elisha Roberson"
# date: "2020-12-08"

# Analysis related to - paley_neurosarcoid_2020_08

# packages
library( here )
library( tidyverse )

# source
source( file = here( "src", "shared_r_functions.R" ) )

# decide graphics parameter based on OS
graphics_device_type <- case_when(
  .Platform$OS.type == "windows" ~ "windows",
  TRUE ~ "cairo"
)

image_width <- 1280
image_height <- 720

# read project info
project_info <- read_tsv( file = here( "info", "unified_project_info.txt" ) )

tcr_summary <- tibble()
tcr_frame <- data.frame()

# process all tcr loop
for ( row_idx in 1:nrow( project_info ) ) {
  rgid <- project_info %>%
    pull( rgid ) %>%
    .[ row_idx ]

  input_tcr <- read_csv( file = here( 'data', 'tcr', rgid, 'filtered_contig_annotations.csv') ) %>%
    filter( full_length == TRUE ) %>%
    filter( productive == TRUE ) %>%
    filter( chain %in% c( 'TRA', 'TRB' ) ) %>%
    mutate( barcode = str_replace( string = barcode, pattern = "\\-1$", replacement = "" ) ) %>%
    mutate( contig_id = str_replace( string = contig_id, pattern = "\\-1", replacement =  "" ) )

  tr_summary <- dplyr::group_by( input_tcr, barcode, chain ) %>%
    dplyr::summarise( count = length( chain ) ) %>%
    pivot_wider( id_cols = barcode, values_from = count, names_from = chain, values_fill = 0 ) %>%
    select( barcode, TRA, TRB )

  tr_table <- with( tr_summary, table( TRA, TRB ) ) %>%
    as.data.frame( . )

  total <- sum( tr_table$Freq )

  tr_table <- mutate( tr_table, Percent = round( ( Freq / total ) * 100.0, 1 ) )

  # VDJ freq
  ggplot_tr_freq <- ggplot( data = tr_table, mapping = aes( TRA, TRB ) ) +
    theme_bw() +
    geom_tile( aes( fill = Freq ), colour = 'white' ) +
    scale_fill_gradient( low = "white", high = "red" ) +
    geom_text( aes( label = Freq ), size = 8 ) +
    gg_bigger_texts

  jpeg( width = image_width,
       height = image_height,
       type = graphics_device_type,
       filename = here( 'results', 'figures', 'tcr_processing', paste0( rgid, "_tcr_freq.jpeg" ) ) )

  print( ggplot_tr_freq )

  dev.off()

  # VDJ percent
  ggplot_tr_percent <- ggplot( data = tr_table, mapping = aes( TRA, TRB ) ) +
    theme_bw() +
    geom_tile( aes( fill = Percent ), colour = 'white' ) +
    scale_fill_gradient( low = "white", high = "red" ) +
    geom_text( aes( label = Percent ), size = 8 ) +
    gg_bigger_texts

  jpeg( width = image_width,
        height = image_height,
        type = graphics_device_type,
        filename = here( 'results', 'figures', 'tcr_processing', paste0( rgid, "_tcr_percent.jpeg" ) ) )

  print( ggplot_tr_percent )

  dev.off()

  # Add summary of counts to big object
  tr_table <- mutate( tr_table, RGID = rgid )

  tcr_summary <- rbind( tcr_summary, tr_table )

  # extract barcodes of good tcr & bad tcr
  tr_good <- tr_summary %>%
    filter( TRB == 1 ) %>%
    filter( TRA > 0 ) %>%
    filter( TRA < 3 ) %>%
    pull( barcode )

  tr_bad <- input_tcr %>%
    pull( barcode ) %>%
    setdiff( x = ., y = tr_good )

  # process good tcr
  number_trs <- select( input_tcr, barcode, chain, contig_id ) %>%
    dplyr::group_by( barcode, chain ) %>%
    dplyr::summarise( chain_number = 1:length( chain ),
                      contig_id = contig_id ) %>%
    as.data.frame( . ) %>%
    select( contig_id, chain_number )

  processed_tcr <- select( input_tcr, contig_id, barcode, chain, raw_clonotype_id, raw_consensus_id, v_gene, d_gene, j_gene, cdr3, cdr3_nt ) %>%
    filter( barcode %in% tr_good ) %>%
    merge( x = number_trs, y = ., by = 'contig_id' ) %>%
    pivot_longer( cols = c( 'raw_consensus_id', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt' ) ) %>%
    mutate( name = paste0( chain, chain_number, "_", name ) ) %>%
    select( -contig_id, -chain, -chain_number ) %>%
    pivot_wider( names_from = name, values_from = value ) %>%
    merge( x = tr_summary, y = ., by = 'barcode' )

  # write qc pass / qc fail
  data.frame( barcode = tr_good, QC = "Pass" ) %>%
    rbind( ., data.frame( barcode = tr_bad, QC = "Fail" ) ) %>%
    write_tsv( x = ., file = here( 'results', 'tcr_qc', paste0( rgid, '_tcr_qc.tsv.gz' ) ) )

  # write processed tcr
  processed_tcr %>%
    write_tsv( x = ., file = here( 'results', 'tcr_meta_data', paste0( rgid, '_tcr_qc.tsv.gz' ) ) )

  tcr_frame <- mutate( processed_tcr, RGID = rgid ) %>%
    dplyr::relocate( RGID, .before = barcode ) %>%
    rbind( tcr_frame, . )
}

# write count summaries
write_tsv( x = tcr_summary, file = here( 'results', 'global_tcr_tr_count.tsv.gz' ) )

# write global tcr frame
write_tsv( x = tcr_frame, file = here( 'results', 'global_tcr_processed.tsv.gz' ) )


# session info
Sys.time()

getwd()

sessionInfo()
