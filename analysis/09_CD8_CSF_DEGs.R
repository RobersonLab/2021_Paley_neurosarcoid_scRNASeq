# title: CD8 CSF diff exp genes
# author: Michael Paley

# DoMultiBarHeatmap source: https://github.com/satijalab/seurat/issues/2201

#load libaries
library( here )
library( ggplot2 )
library( dplyr )
library( tidyr )
library( readr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
library( sctransform )
library( RColorBrewer )
library( grid )
library( escape )

# Load the Seurat Object 'scd' ("Single Cell Data")
scd <- readRDS( file = here( 'results', 'objects', 'scd.rds' ) )

scd # Sanity check
head( Idents(object = scd) ) # Sanity check

#Load Metadata
MetaData <- read_tsv( file = here( 'info', 'unified_project_info.tsv') )
nfiles <- nrow( MetaData )

# source (Important for DoMultiBarHeatmap)
source( file = here( 'src', 'shared_r_functions.R' ) )

# Color scale for heatmaps
mapal <- colorRampPalette( RColorBrewer::brewer.pal( 11, "RdBu" ) )( 256 )

### MISC SECTION: ADD FIGURE OF HIGHLIGHTED MAIN LINEAGES ###
Idents(scd) <- 'lin_1'
Idents(scd) <- factor(Idents(scd), levels =  sort( levels( Idents( scd ) ) ) ) # Needed to alphabetize the cluster level
highlight <- list() # create a list of the lineages to highlight
for ( i in c( 'B' , 'CD4T' , 'CD8T' , 'M' , 'NK' ) ) {
    highlight[[i]] <- WhichCells( scd, idents = i )
}
lin.colors <- c( '#D682FD' , '#FF9200' , '#7A80FD' , '#FF7E78' , '#008F91' ) # This matches the alphabetical order of the highlighted lineages
umap.plot <- DimPlot(scd, group.by = 'lin_1',
                     cells.highlight = highlight,
                     cols.highlight = rev(lin.colors),
                     sizes.highlight = 0.1,
                     pt.size = 0.025) + # application of colors between DimPlot and DoHeatmap is opposite
        ggtitle( 'Main Lineages' )
file.name <- 'UMAP.annotated.main.lin.highlighted.pdf'

pdf( file = here( 'results', 'figures', file.name ) )
print( umap.plot )
dev.off()
### END OF MISC SECTION ###

# Three comparisons for differential gene expression
# 1) Blood CD8s vs CSF CD8s in NS subjects (Non-Naive CD8s)
### a) with TRM cluster
### b) without TRM cluster
# 2) NS vs control CD8s in the CSF
### a) with TRM cluster
### b) without TRM cluster
# 3) NS001 Active vs Remission in the CSF

#########################################
# Annotate / Subset CD8 T cell clusters #
#########################################
# Read in Lineage Annotation File and Add Annotation as MetaData to Seurat Object
Idents(scd) <- 'integrated_snn_res.2.5'
lin_annotation <- read_tsv( file = here( 'info', 'lin_annotation2.tsv' ) )
clusters.2.5 <- as.data.frame( Idents(scd) )
colnames(clusters.2.5) <- "cluster"
clusters.2.5$cluster <- as.numeric(as.character(clusters.2.5$cluster)) # This step was needed because the Idents were Factors that had Levels that did not correspond numerically
new.cluster.ids <- left_join( clusters.2.5, lin_annotation, by = "cluster" )
rownames(new.cluster.ids) <- rownames(clusters.2.5)
scd <- AddMetaData( scd, new.cluster.ids$lineage2, col.name = 'lin_2')

subjects <- unlist( strsplit( scd$Sample , "_" ) )[ c(TRUE,FALSE) ] %>%
            as.data.frame()
colnames(subjects) <- c( 'Subjects' )
rownames(subjects) <- rownames(clusters.2.5)
scd <- AddMetaData( scd, subjects, col.name = 'Subjects')

# Print UMAP of 'lin_2' annotations highlighted
section <- "CD8.T.Clusters"
Idents(scd) <- 'lin_2'
Idents(scd) <- factor(Idents(scd), levels =  sort( levels( Idents( scd ) ) ) ) # Needed to re-sort the cluster level
highlight <- list() # create a list of the CD8T clusters to highlight
for ( i in c( 1:10 ) ) {
    j <- (i - 1)
    cd8.cluster <- paste0( 'CD8T_', j )
    highlight[[cd8.cluster]] <- WhichCells( scd, idents = cd8.cluster )
}

cd8.cluster.colors <- brewer.pal(10, 'Paired')
umap.plot <- DimPlot(scd,
                     group.by = 'lin_2',
                     cells.highlight = highlight,
                     cols.highlight = rev(cd8.cluster.colors),
                     sizes.highlight = 0.1,
                     pt.size = 0.025) + # application of colors between DimPlot and DoHeatmap is opposite
        ggtitle( section )
file.name <- paste0(section, '.pdf')

pdf( file = here( 'results', 'figures', 'CD8T', file.name ) )
print( umap.plot )
dev.off()

# Subset CD8 T cells and generate markers for each cluster
Idents(scd) <- 'lin_1'
tsub <- subset( scd, idents = 'CD8T' ) # Subset CD8 T cells

Idents(tsub) <- 'lin_2' # Label with CD8T clusters
Idents(tsub) <- factor(Idents(tsub), levels =  sort( levels( Idents( tsub ) ) ) ) # Needed to re-sort the cluster levels

# NORMALiZE AND SCALL DATA ON SUBSET
DefaultAssay(tsub) <- "RNA"
tsub <- NormalizeData(tsub) # Normalize data
all.genes <- rownames(tsub)
tsub <- ScaleData(tsub, features = all.genes) # Scale data
cd8.markers <- FindAllMarkers(tsub,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

x <- filter( cd8.markers, p_val_adj < 0.05 ) %>%
    dplyr::group_by( cluster ) %>%
    arrange( desc( avg_log2FC ), .by_group = TRUE )

file.name <- paste0( section, '.tsv' )
write.table( x = x,
             file = here( 'results', 'DEGs', file.name ),
             sep = "\t",
             quote = FALSE,
             row.names = FALSE )

# Heatmap
Idents(tsub) <- 'Subjects'
subject.colors <- brewer.pal(8, 'Set2')
colors.to.use <- list( lin_2 = cd8.cluster.colors,  Tissue = c( 'blue' , 'red' ), Subjects = subject.colors )

cd8.markers.heatmap <- DoMultiBarHeatmap(tsub,
                                         features = x$gene,
                                         group.by = 'lin_2',
                                         additional.group.by = c('Tissue', 'Subjects'),
                                         additional.group.sort.by = c('Tissue', 'Subjects'),
                                         group.bar = TRUE,
                                         cols.use = colors.to.use,
                                         assay = 'RNA',
                                         label = TRUE,
                                         size = 5,
                                         draw.lines = TRUE,
                                         lines.width = 20) +
    NoLegend() +
    theme( axis.text.y = element_text(size = 0.5) ) + # decrease size of gene names
    scale_fill_gradientn( colours = rev(mapal) ) + # change color gradient to Red-Blue
    theme( plot.margin = unit( c(2,2,2,2) ,"cm") )

file.name <- paste0(section, '.heatmap.pdf')
pdf( file = here( 'results', 'figures', 'CD8T', file.name ),
     width = 8,
     height = 11 )
print( cd8.markers.heatmap )
dev.off()

rm( cd8.markers.heatmap, x )

##################################################################
# 1 Find markers CSF CD8 T cells in NS subjects (Non-Naive CD8s) #
##################################################################
section <- "Non_Naive_CD8_NS_blood_vs_csf"
a <- 'with_TRM'
Idents(tsub) <- 'lin_2'
non.naive.cd8s <- paste0( 'CD8T_', 1:9 ) # Manual list of Non-Naive CD8 T Cell Clusters
non.naive.sub <- subset( tsub, idents = non.naive.cd8s ) # Subset Non-naive CD8 T cells

Idents(non.naive.sub) <- 'Status' # Change Idents for active, inactive, & control
non.naive.sub <- subset( non.naive.sub, idents = 'active' ) # Subset CD8 T cells for active NS, e.g. removes cells from inactive and controls

Idents(non.naive.sub) <- 'Tissue'
DefaultAssay(non.naive.sub) <- "RNA"
tissue.markers <- FindAllMarkers(non.naive.sub,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)

x <- filter( tissue.markers, p_val_adj < 0.05 ) %>%
    group_by(cluster)
file.name <- paste( section, a, 'tsv', sep = '.' )
write.table( x = x,
             file = here( 'results', 'DEGs', file.name ),
             sep = "\t",
             quote = FALSE,
             row.names = FALSE )

#Heatmap
non.naive.cluster.colors <- cd8.cluster.colors[2:10]
colors.to.use <- list( lin_2 = non.naive.cluster.colors,
                       Tissue = c( 'blue' , 'red' ),
                       Subjects = subject.colors[1:6] ) # change number of subjects to only NS

cd8.markers.heatmap <- DoMultiBarHeatmap(non.naive.sub,
                                         features = x$gene,
                                         group.by = 'Tissue',
                                         additional.group.by = c('lin_2', 'Subjects'),
                                         additional.group.sort.by = c('lin_2', 'Subjects'),
                                         group.bar = TRUE,
                                         cols.use = colors.to.use,
                                         assay = 'RNA',
                                         label = TRUE,
                                         size = 5,
                                         draw.lines = TRUE,
                                         lines.width = 20) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
    scale_fill_gradientn(colours = rev(mapal)) + # change color gradient to Red-Blue
    theme( plot.margin = unit( c( 25, 75, 10, 10 ) ,"pt") ) # margin order is top, right, bottom, left

file.name <- paste0(section, '.', a, '.heatmap.pdf')
pdf( file = here( 'results', 'figures', 'CD8T', file.name ),
     width = 8,
     height = 11 )
print( cd8.markers.heatmap )
dev.off()

rm( non.naive.sub, tissue.markers, x, cd8.markers.heatmap )

# Repeat without the 'TRM' / 'CD8T_9' cluster
b <- 'no_TRM'
non.naive.cd8s <- paste0( 'CD8T_', 1:8 ) # Manual list of Non-Naive CD8 T Cell Clusters
non.naive.sub <- subset( tsub, idents = non.naive.cd8s ) # Subset Non-naive CD8 T cells
Idents(non.naive.sub) <- 'Tissue'
DefaultAssay(non.naive.sub) <- "RNA"
markers <- FindAllMarkers(non.naive.sub,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

x <- filter( markers, p_val_adj < 0.05 ) %>%
    group_by(cluster)

file.name <- paste( section, b, 'tsv', sep = '.' )
write.table( x = x,
            file = here( 'results', 'DEGs', file.name ),
            sep = '\t',
            quote = FALSE,
            row.names = FALSE )

#Heatmap
non.naive.cluster.colors <- cd8.cluster.colors[2:9]
colors.to.use <- list( lin_2 = non.naive.cluster.colors,
                       Tissue = c( 'blue' , 'red' ),
                       Subjects = subject.colors[1:6] ) # change number of subjects to only NS
cd8.markers.heatmap <- DoMultiBarHeatmap(non.naive.sub,
                                         features = x$gene,
                                         group.by = 'Tissue',
                                         additional.group.by = c('lin_2', 'Subjects'),
                                         additional.group.sort.by = c('lin_2', 'Subjects'),
                                         group.bar = TRUE,
                                         cols.use = colors.to.use,
                                         assay = 'RNA',
                                         label = TRUE,
                                         size = 5,
                                         draw.lines = TRUE,
                                         lines.width = 20) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
    scale_fill_gradientn(colours = rev(mapal)) + # change color gradient to Red-Blue
    theme( plot.margin = unit( c( 25, 75, 10, 10 ) ,"pt") ) # margin order is top, right, bottom, left

file.name <- paste0(section, '.', b, '.heatmap.pdf')
pdf( file = here( 'results', 'figures', 'CD8T', file.name ),
     width = 8,
     height = 11 )
print( cd8.markers.heatmap )
dev.off()

rm( non.naive.sub, markers, x )

############################################
# 2 Find markers for NS in CSF CD8 T cells #
############################################
# Identify CD8 T cell markers specific to NS; NS vs control CD8s in the CSF
section <- "CSF_CD8_NS_vs_ctrl"
a <- 'with_TRM'
Idents(tsub) <- 'Tissue'
csf.sub <- subset( tsub, idents = 'CSF' ) # Subset CSF CD8 T cells
Idents(csf.sub) <- 'Status' # Change Idents for active, inactive, & control
csf.sub <- subset(csf.sub, idents = c('active', 'control') ) # Subset CD8 T cells in active NS and controls, e.g. removes cells from inactive
DefaultAssay(csf.sub) <- "RNA"
markers <- FindAllMarkers(csf.sub,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

x <- filter( markers, p_val_adj < 0.05 ) %>%
    group_by(cluster)

file.name <- paste( section, a, 'tsv', sep = '.' )
write.table( x = x,
            file = here( 'results', 'DEGs', file.name ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE )

#Heatmap
colors.to.use <- list( lin_2 = cd8.cluster.colors,
                       Status = c( 'red', 'gray' ),
                       Subjects = subject.colors[c(1,3:8)] ) # select colors not used by 'inactive' subject
cd8.markers.heatmap <- DoMultiBarHeatmap(csf.sub,
                                         features = x$gene,
                                         group.by = 'Status',
                                         additional.group.by = c('lin_2', 'Subjects'),
                                         additional.group.sort.by = c('lin_2', 'Subjects'),
                                         group.bar = TRUE,
                                         cols.use = colors.to.use,
                                         assay = 'RNA',
                                         label = TRUE,
                                         size = 5,
                                         draw.lines = TRUE,
                                         lines.width = 20) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
    scale_fill_gradientn(colours = rev(mapal)) + # change color gradient to Red-Blue
    theme( plot.margin = unit( c( 25, 75, 10, 10 ) ,"pt") ) # margin order is top, right, bottom, left

file.name <- paste0(section, '.', a, '.heatmap.pdf')
pdf( file = here( 'results', 'figures', 'CD8T', file.name ),
     width = 8,
     height = 11 )
print( cd8.markers.heatmap )
dev.off()

rm( markers, x )

# Repeat without the 'TRM' / 'CD8T_9' cluster
b <- 'no_TRM'

Idents(csf.sub) <- 'lin_2'
non.trm.cd8s <- paste0( 'CD8T_', 0:8 ) # Manual list of Non-TRM CD8 T Cell Clusters
non.trm.csf.sub <- subset( csf.sub, idents = non.trm.cd8s ) # Subset Non-TRM CD8 T cells
Idents(non.trm.csf.sub) <- 'Status'
DefaultAssay(non.trm.csf.sub) <- "RNA"
markers <- FindAllMarkers(non.trm.csf.sub,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

x <- filter( markers, p_val_adj < 0.05 ) %>%
    group_by(cluster)

file.name <- paste( section, b, 'tsv', sep = '.' )
write.table( x = x,
             file = here( 'results', 'DEGs', file.name ),
             sep = "\t",
             quote = FALSE,
             row.names = FALSE )

#Heatmap
colors.to.use <- list( lin_2 = cd8.cluster.colors[1:9],
                       Status = c( 'red', 'gray' ),
                       Subjects = subject.colors[c(1,3:8)] ) # select colors not used by 'inactive' subject and CD8T clusters

cd8.markers.heatmap <- DoMultiBarHeatmap(non.trm.csf.sub,
                                         features = x$gene,
                                         group.by = 'Status',
                                         additional.group.by = c('lin_2', 'Subjects'),
                                         additional.group.sort.by = c('lin_2', 'Subjects'),
                                         group.bar = TRUE,
                                         cols.use = colors.to.use,
                                         assay = 'RNA',
                                         label = TRUE,
                                         size = 5,
                                         draw.lines = TRUE,
                                         lines.width = 20 ) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
    scale_fill_gradientn(colours = rev(mapal)) + # change color gradient to Red-Blue
    theme( plot.margin = unit( c( 25, 75, 10, 10 ) ,"pt") ) # margin order is top, right, bottom, left

file.name <- paste0(section, '.', b, '.heatmap.pdf')
pdf( file = here( 'results', 'figures', 'CD8T', file.name ),
     width = 8,
     height = 11 )
print( cd8.markers.heatmap )
dev.off()

rm( non.trm.csf.sub, markers, x )

###############################
# 3 NS001 Active vs Remission #
###############################
section <- "NS001_CSF_CD8_active_vs_inactive"
Idents(tsub) <- 'Tissue'
csf.sub <- subset( tsub, idents = 'CSF' ) # Subset CSF CD8 T cells
Idents(csf.sub) <- 'Subjects'
ns001.sub <- subset( csf.sub, idents = c( 'NS001-1', 'NS001-2') ) # Subset CSF CD8 T cells from NS001 before and after treatment
DefaultAssay(ns001.sub) <- "RNA"
markers <- FindAllMarkers(ns001.sub,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

x <- filter( markers, p_val_adj < 0.05 ) %>%
    group_by(cluster)

file.name <- paste( section, 'tsv', sep = '.' )
write.table( x = x,
            file = here( 'results', 'DEGs', file.name ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE )

#Heatmap
colors.to.use <- list( lin_2 = cd8.cluster.colors,
                       Status = c( 'red', 'gray' ),
                       Subjects = subject.colors[1:2] ) # select colors used by subject
cd8.markers.heatmap <- DoMultiBarHeatmap(ns001.sub,
                                         features = x$gene,
                                         group.by = 'Status',
                                         additional.group.by = c('lin_2', 'Subjects'),
                                         additional.group.sort.by = c('lin_2', 'Subjects'),
                                         group.bar = TRUE, cols.use = colors.to.use,
                                         assay = 'RNA',
                                         label = TRUE,
                                         size = 5,
                                         draw.lines = TRUE,
                                         lines.width = 20) +
    NoLegend() +
    theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
    scale_fill_gradientn(colours = rev(mapal)) + # change color gradient to Red-Blue
    theme( plot.margin = unit( c( 25, 75, 10, 10 ) ,"pt") ) # margin order is top, right, bottom, left

file.name <- paste0(section, '.heatmap.pdf')
pdf( file = here( 'results', 'figures', 'CD8T', file.name ),
     width = 8,
     height = 11 )
print( cd8.markers.heatmap )
dev.off()

########################################################################
# Calculating IFN score for CD8 T cell clusters per Tissue per Subject #
########################################################################
# Escape - ssGSEA
# source: https://ncborcherding.github.io/vignettes/escape_vignette.html#1_loading_libraries
if ( dir.exists( here( 'results', 'Gene.Sets' ) ) ) {
    message( "Gene.Sets exists and is a directory.\n" )
} else {
    message( "Gene.Sets does not exist - creating.\n" )
    dir.create( path = here( 'results', 'Gene.Sets' ) )
}

# Create a 'Tissue-Diagnosis' column of metadata
tsub@meta.data$Tissue.GeneralStatus <- paste(tsub@meta.data$Tissue,
                                             tsub@meta.data$GeneralStatus,
                                             sep = '.')

# Get gene sets
GS <- getGeneSets(library = "H") # Hallmark

# Calculate GSEA for tsub; choose n for groups based on number of cells in smallest sample
ES <- enrichIt(obj = tsub,
               gene.sets = GS,
               groups = 1000,
               cores = 4) # Groups are the number of cells to separate the enrichment calculation; cores are for parallelization

# Add Enrichment scores to tsub
tsub <- AddMetaData(tsub, ES)

#Calculate Statistical Significance
ES2 <- data.frame( tsub[[]][ , (ncol(tsub[[]]) - ncol(ES) + 1):ncol(tsub[[]]) ],
                   GeneralStatus = tsub$GeneralStatus,
                   Tissue = tsub$Tissue,
                   Tissue.GeneralStatus = tsub$Tissue.GeneralStatus )

output.anova <- getSignificance(ES2,
                                group = "Tissue.GeneralStatus",
                                fit = "ANOVA")  # fit =  linear.model, ANOVA, or T.test ; linear.model is based on limma

x <- data.frame(Gene_Set = rownames(output.anova), output.anova)

file.name <- paste( 'CD8T.ssGSEA.Stats.ANOVA', 'tsv', sep = '.' )
write.table( x = x,
             file = here( 'results', 'Gene.Sets', file.name ),
             sep = "\t",
             quote = FALSE,
             row.names = FALSE )

# Visualize Gene Sets with Violin Plots custom y scale
plot_list <- list()
for ( i in 1:length(GS) ) {
    plot_list[[i]] <- VlnPlot( object = tsub,
                               features = GS[[i]]@setName,
                               cols = c( '#939598' , '#FF4B20' ),
                               pt.size = 0,
                               assay = 'RNA',
                               split.by = 'GeneralStatus' ) +
        theme( text = element_text(size = 16),
               axis.text.x = element_text(angle = 0, hjust = 0.5) ) +
        ylim( 0, 0.8 )
}

for ( i in 1:length(GS) ) {
    file.name <- paste( 'Reg', 'Vln', 'Plot', 'CD8T', GS[[i]]@setName, 'pdf', sep = '.' )
    pdf( file = here( 'results', 'Gene.Sets', file.name ),
         width = 8,
         height = 11 )
    print( plot_list[[i]] )
    dev.off()
}

# Save the objects
saveRDS( scd, file = here( 'results', 'objects', 'scd.rds' ) )

# session info
Sys.time()
getwd()
sessionInfo()
