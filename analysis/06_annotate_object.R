# title: annotate object
# author: Michael Paley

# load libaries
library( here )
library( ggplot2 )
library( dplyr )
library( tidyr )
library( tibble )
library( Seurat )
library( sctransform )

# Load the Seurat Object 'scd' ("Single Cell Data")
scd <- readRDS( file = here( 'results', 'objects', 'scd.rds' ) )

# sanity checks
scd
head( Idents( object = scd ) )

# Read in Lineage Annotation File and Add Annotation as MetaData to Seurat Object then convert to Idents
lin_annotation <- read.delim( file = here( 'info', 'lin_annotation1.tsv' ) )

chosen_clusters <- as.data.frame( Idents( scd ) )

colnames( chosen_clusters ) <- "cluster"
chosen_clusters$cluster <- as.integer( as.character( chosen_clusters$cluster ) ) # This step was needed because the Idents were Factors that had Levels that did not correspond numerically

new.cluster.ids <- left_join( chosen_clusters, lin_annotation, by = "cluster" )
rownames( new.cluster.ids ) <- rownames( chosen_clusters )
scd <- AddMetaData( scd, new.cluster.ids$lineage, col.name = 'lin_1' )
Idents( scd ) <- "lin_1"

# Print annotated UMAP plot
pdf( file = here( 'results', 'figures', paste( 'UMAP', 'annotated', 'pdf', sep = '.' ) ) )
print(
  DimPlot( scd,
           group.by = 'ident',
           pt.size = 0.1)
)
dev.off()

# Save the object
saveRDS( object = scd, file = here( 'results', 'objects', 'scd.rds' ) )

# session info
Sys.time()
getwd()
sessionInfo()
