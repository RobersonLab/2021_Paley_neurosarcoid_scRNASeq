# title: UMAP and find clusters
# author: Michael Paley

# load libaries
library( here )
library( ggplot2 )
library( cowplot )
library( dplyr )
library( tidyr )
library( tibble )
library( Seurat )
library( sctransform )
library( patchwork )

# Load the Seurat Object 'scd' ("Single Cell Data")
load( file = here( 'results', 'objects', 'integrated_data.RData' ) )
scd <- integrated_data
rm( integrated_data )

# Remove TCR and BCR genes from "Variable Features" in Seurat Object
tcr_bcr_gene_grep_string <- "^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC|^IGHV|^IGHD|^IGHJ|^IGHM|^IGHE|^IGHG|^IGHA|^IGLV|^IGLD|^IGLC|^IGKV|^IGKD|^IGKC"
VariableFeatures( scd ) <- VariableFeatures( scd )[!grepl(tcr_bcr_gene_grep_string, scd@assays$integrated@var.features)]

# Perform PCA
scd <- RunPCA( scd, verbose = FALSE )
message( "RunPCA complete" )

# Perform Dimensionality Reduction with UMAP
scd <- RunUMAP( scd, dims = 1:30, verbose = FALSE )
message( "RunUMAP complete" )

scd <- FindNeighbors( scd, dims = 1:30, verbose = FALSE )
message( "FindNeighbors complete" )

# Find Clusters based on a defined resolution (RES)
RES <- 2.5
scd <- FindClusters( scd, resolution = RES, verbose = FALSE )
message( "FindClusters complete" )

# Save the object
saveRDS( object = scd, file = here( 'results', 'objects', 'scd.rds' ) )
message( "Save Object complete" )

# Print Plots grouped by MetaData
n <- colnames( scd@meta.data[,4:8] )
plot_list <- list()
for ( i in 1:length( n ) ) {
    plot_list[[i]] <- DimPlot( scd,
                               group.by = n[i],
                               pt.size = 0.025 )
}

for ( i in 1:length(n) ) {
    pdf( file = here( 'results', 'figures', paste( 'UMAP', n[i], 'pdf', sep = '.' ) ) )
    print( plot_list[[i]] )
    dev.off()
}
message( "Printed UMAP grouped by MetaData" )

# Print Plots of Lineage Markers #
DefaultAssay( scd ) <- "SCT"
lin.markers <- read.delim( file = here( 'info', 'lin_markers.tsv' ) )

n <- colnames(lin.markers)
plot_list <- list()
for ( i in 1:length(n) ) {
    plot_list[[i]] <- FeaturePlot( scd,
                                   features = lin.markers[[i]],
                                   pt.size = 0.1,
                                   ncol = 2 )
}

for ( i in 1:length(n) ) {
    pdf( file = here( 'results', 'figures', paste( 'UMAP', 'lineage', n[i], 'pdf', sep = '.' ) ) )
    print( plot_list[[i]] )
    dev.off()
}
message( "Printed UMAP based on Lineage Markers" )

# session info
Sys.time()
getwd()
sessionInfo()
