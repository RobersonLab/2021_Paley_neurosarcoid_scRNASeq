# title: plot cluster frequency
# author: Michael Paley

# load libaries
library( here )
library( ggplot2 )
library( dplyr )
library( tidyr )
library( tibble )
library( stringr )
library( Seurat )
library( sctransform )
library( RColorBrewer )

# source it
source( file = here( 'src', 'shared_r_functions.R' ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
scd <- readRDS( file = here( 'results', 'objects', 'scd.rds' ) )

scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Set color scheme base for plots
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]

col_vector <- unlist(
  mapply(
    brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
    )
  )

# Set vectors of Samples, Subjects, and Tissue
Idents(scd) <- 'Sample'
nsamples <- levels( Idents(scd) )
nsplit <- str_split(nsamples, '_', simplify = TRUE) %>%
  as.data.frame()
colnames(nsplit) <- c('Subject', 'Tissue')
nsplit <- arrange(nsplit, Tissue)
njoin <- paste0(nsplit$Subject,'_',nsplit$Tissue)

# Subset CD8T cells
section <- 'Tcell.CD8'
Idents(scd) <- 'lin_1'
tsub <- subset(scd, idents = c('CD8T') )

# Calculate frequencies of T cell clusters
Idents(tsub) <- 'integrated_snn_res.2.5'
n <- length( levels( Idents(tsub) ) )
counts <- table( tsub$Sample, Idents(tsub) )
freq <- round( counts / rowSums( counts ), digits = 4 ) #normalize to percentages within each cell state

# Save calulations
write.csv( x = as.matrix( counts ),
           file = here( 'results', 'calculations', paste0( section, '.counts.per.sample.matrix.csv' ) ),
           quote = FALSE,
           row.names = TRUE )

write.csv( x = as.matrix( freq ),
           file = here( 'results', 'calculations', paste0( section,'.freq.per.sample.matrix.csv' ) ),
           quote = FALSE,
           row.names = TRUE )

# Convert Cluster Frequencies Table to Dataframe for plotting
counts <- as.data.frame( counts )
colnames( counts ) <- c( 'Sample', 'Cluster', 'Counts' )

freq <- as.data.frame(freq)
colnames( freq ) <- c( 'Sample', 'Cluster', 'Freq' )

# Save calulations
write.csv( x = as.data.frame(counts),
           file = here( 'results', 'calculations', paste0( section,'.counts.per.sample.df.csv' ) ),
           quote = FALSE,
           row.names = TRUE )

write.csv( x = as.data.frame(freq),
           file = here( 'results', 'calculations', paste0( section,'.freq.per.sample.df.csv' ) ),
           quote = FALSE,
           row.names = TRUE )

# Print Cluster Frequencies in Stacked Bar Plot
pdf( file = here(
  'results',
  'calculations',
  paste( section, 'freq', 'per', 'sample',' pdf', sep = '.' ) ) )

plot_out <- ggplot(data = counts,
                   aes(
                     fill=Cluster,
                     y=Counts,
                     x=factor( Sample, level = njoin ) ) ) +
  geom_bar(position="fill", stat="identity") +
  xlab( "Sample" ) +
  ylab( "Proportion" ) +
  scale_fill_manual( values = sample( col_vector, n ), aesthetics = "fill" ) + # Using 'n' defined above
  theme(
    panel.background = element_rect( fill="white" ),
    axis.text.x = element_text(
      angle=45,
      hjust = 1,
      vjust=1,
      size = 12,
      face = "bold" ),
    axis.text.y = element_text(
      size = 12,
      face = "bold" ) )

print( plot_out )

dev.off()

message( paste0( section, " ", "section complete" ) )

# session info
Sys.time()
getwd()
sessionInfo()
