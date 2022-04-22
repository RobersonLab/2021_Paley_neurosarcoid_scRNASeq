# title: looking at clones across tissues
# author: Michael Paley

# Calculate Average Expression and create Rule to automatically ID Lineages

#load libaries
library( here )
library( tidyverse )
#library( ggplot2 )
library( cowplot )
#library( dplyr )
#library( tidyr )
#library( purrr )
#library( tibble )
#library( stringr )
library( Seurat )
library( pheatmap )
library( RColorBrewer )
library( SeuratWrappers )
library( monocle )
library( VGAM )

# Load the Seurat Object 'scd' ("Single Cell Data")
scd <- readRDS( file = here( 'results', 'objects', 'scd.rds' ) )
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

#Load Metadata
MetaData <- read_tsv( file = here( 'info', 'unified_project_info.tsv' ) ) %>%
  mutate( Sample = rgsm ) %>%
  mutate( Subject = individual )
nfiles <- nrow( MetaData )

# Create TCR clonotypes via TRBV-CDR3 & TRAV-CDR3
# unique.clones <- data.frame( barcode = colnames(scd),
#                              Sample = scd$Sample,
#                              clone_lin_2 = scd$clone_lin_2,
#                              Tissue = scd$Tissue )

unique.clones <- data.frame( barcode = colnames(scd),
                             Sample = scd$Sample,
                             clone_lin_2 = scd$lin_2,
                             Tissue = scd$Tissue )

unique.clones$tcr.seq <- paste( scd$TRB1_v_gene,
                                scd$TRB1_cdr3,
                                scd$TRA1_v_gene,
                                scd$TRA1_cdr3, sep = '_' )

unique.clones <- unique.clones %>%
  filter( !is.na(clone_lin_2) )

# unique.clones <- left_join(unique.clones,
#                            MetaData[,c('Sample','Subject')],
#                            by = 'Sample')

unique.clones <- left_join( unique.clones,
                            select( MetaData, Sample, Subject ),
                            by = 'Sample')

# cd4.clones <- unique.clones %>%
#   filter(clone_lin_2 == "CD4")
#
# cd8.clones <- unique.clones %>%
#   filter(clone_lin_2 == "CD8")

cd4.clones <- unique.clones %>%
  filter( str_detect( string = clone_lin_2, pattern = "^CD4" ) )

cd8.clones <- unique.clones %>%
  filter( str_detect( string = clone_lin_2, pattern = "^CD8" ) )

# Visualize clonotypes in paired samples for CD4s and CD8s
# For each clone, (1) subset based on patient/subject, (2) table clone and tissue, (3) calculate % of each clonotype, (4) plot x vs y of freq in each tissue

# CD4 T CELL LOOP #
cd4.enriched <- c()
cd4.clones.enriched <- list()
for (i in unique(MetaData$Subject) ) {
  pt.clones <- subset(cd4.clones, Subject == i)

  pt.clone.table <- table(pt.clones$tcr.seq, pt.clones$Tissue) %>%
    as.data.frame()

  colnames(pt.clone.table) <- c('tcr.seq', 'Tissue', 'Count')
  pt.clone.wide <- pivot_wider( data = pt.clone.table,
                                names_from = Tissue,
                                values_from = Count) %>%
    as.data.frame()

  pt.clone.pct <- mutate(pt.clone.wide,
                         CSF_pct = 100 * CSF / sum(CSF),
                         PBMC_pct = 100 * PBMC / sum(PBMC) )

  # Plot CSF-enriched (>5x) and expanded clonotypes (n > 5)
  pt.clone.pct$CSF_PBMC_ratio <- pt.clone.pct$CSF_pct / pt.clone.pct$PBMC_pct
  pt.clone.enriched <- pt.clone.pct %>%
    filter( CSF_PBMC_ratio > 5 ) %>%
    filter( CSF > 5)

  cd4.enriched[i] <- nrow(pt.clone.enriched)
  cd4.clones.enriched[[i]] <- pt.clone.enriched$tcr.seq

  # Remove Intermediate objects
  rm( pt.clones, pt.clone.table, pt.clone.wide, pt.clone.pct, pt.clone.enriched )
}

# CD8 T CELL LOOP #
cd8.enriched <- c()
cd8.clones.enriched <- list()
for (i in unique(MetaData$Subject) ) {
  pt.clones <- subset(cd8.clones, Subject == i)

  pt.clone.table <- table(pt.clones$tcr.seq, pt.clones$Tissue) %>%
    as.data.frame()

  colnames(pt.clone.table) <- c('tcr.seq', 'Tissue', 'Count')

  pt.clone.wide <- pivot_wider( data = pt.clone.table,
                                names_from = Tissue,
                                values_from = Count ) %>%
    as.data.frame()

  pt.clone.pct <- mutate(pt.clone.wide,
                         CSF_pct = 100 * CSF / sum(CSF),
                         PBMC_pct = 100 * PBMC / sum(PBMC) )

  # Plot CSF-enriched (>5x) and expanded clonotypes (n > 5)
  pt.clone.pct$CSF_PBMC_ratio <- pt.clone.pct$CSF_pct / pt.clone.pct$PBMC_pct

  pt.clone.enriched <- pt.clone.pct %>%
    filter( CSF_PBMC_ratio > 5 ) %>%
    filter( CSF > 5 )

  cd8.enriched[i] <- nrow(pt.clone.enriched)
  cd8.clones.enriched[[i]] <- pt.clone.enriched$tcr.seq

  csv.file.name <- paste( 'CD8','clones','CSF','vs','PBMC', i,'csv', sep='.' )
  write.csv( x = pt.clone.pct,
             file = here( 'results', 'calculations', 'VDJ', csv.file.name ) )

  # Remove Intermediate objects
  rm( pt.clones, pt.clone.table, pt.clone.wide, pt.clone.pct, pt.clone.enriched )
}

enriched.cd8.tcrs <- unlist(cd8.clones.enriched,
                            recursive = TRUE,
                            use.names = FALSE) %>%
  as.factor( . ) %>%
  as.character( . )

# Create new column of clonotypes based on tcr.seq data
scd$tcr.seq <- paste(scd$TRB1_v_gene,
                     scd$TRB1_cdr3,
                     scd$TRA1_v_gene,
                     scd$TRA1_cdr3,
                     sep = '_')
Idents(scd) <- 'tcr.seq'

# Save the barcodes assocaited with CSF-enriched TCRs
highlight <- WhichCells(scd, idents = enriched.cd8.tcrs)

enriched.cd8.plot <- DimPlot(scd,
                             group.by = 'Sample',
                             cells.highlight = highlight,
                             cols.highlight = '#FF0000',
                             sizes.highlight = 0.01,
                             pt.size = 0.001) +
  ggtitle( 'CSF-enriched CD8 T Cell Clonotypes' ) +
  theme(legend.position = 'none')

file.name <- paste('CSF-enriched_CD8_T_Cell_Clonotypes', 'pdf', sep='.')
pdf( file = here( 'results', 'figures', 'VDJ', file.name ),
     width = 8,
     height = 8 )
print(enriched.cd8.plot)
dev.off()

# Save TCRs & Barcodes of CSF-enriched TCRs
csv.file.name <- paste('CSF','specific','TCR','seq','csv', sep='.')
write.csv( x = enriched.cd8.tcrs,
           file = here( 'results', 'calculations', 'VDJ', csv.file.name ) )

csv.file.name <- paste('CSF','specific','TCR','barcodes','csv', sep='.')
write.csv( x = highlight,
           file = here( 'results', 'calculations', 'VDJ', csv.file.name ) )

# Subset the data based on the CSF-enriched TCRs
tcr.sub <- subset(scd, subset = tcr.seq %in% enriched.cd8.tcrs)

# Calculate the number of CSF-enriched T cells in each cluster per donor
# Set color scheme base for plots
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Set vectors of Samples, Subjects, and Tissue
section <- 'CSF-enriched.CD8.TCRs'
Idents(tcr.sub) <- 'Sample'
nsamples <- levels( Idents(tcr.sub) )
nsplit <- str_split(nsamples, '_', simplify = TRUE) %>%
  as.data.frame( . )
colnames(nsplit) <- c('Subject', 'Tissue')
nsplit <- arrange(nsplit, Tissue)
njoin <- paste0(nsplit$Subject,'_',nsplit$Tissue)
# Calculate frequencies of clusters per sample for CSF-enriched TCRs
Idents(tcr.sub) <- 'integrated_snn_res.2.5'
n <- length( levels( Idents(tcr.sub) ) )
counts <- table( tcr.sub$Sample, Idents(tcr.sub) )
freq <- round(counts/rowSums(counts), digits = 4) #normalize to percentages within each cell state

# Save calulations
write.csv( x = as.matrix( counts ),
           file = here( 'results',
                        'calculations',
                        'VDJ',
                        paste0(section, '.counts.per.sample.matrix.csv') ),
           quote = FALSE,
           row.names = TRUE )

write.csv( x = as.matrix( freq ),
           file = here( 'results',
                        'calculations',
                        'VDJ',
                        paste0(section,'.freq.per.sample.matrix.csv') ),
           quote = FALSE,
           row.names = TRUE )

# Convert Cluster Frequencies Table to Dataframe for plotting
counts <- as.data.frame(counts)
colnames(counts) <- c('Sample', 'Cluster', 'Counts')
freq <- as.data.frame(freq)
colnames(freq) <- c('Sample', 'Cluster', 'Freq')

# Save calulations
write.csv( as.data.frame(counts),
           file = here( 'results',
                        'calculations',
                        'VDJ',
                        paste0(section,'.counts.per.sample.df.csv' ) ),
           quote = FALSE,
           row.names = TRUE )

write.csv( x = as.data.frame( freq ),
           file = here( 'results',
                        'calculations',
                        'VDJ',
                        paste0(section,'.freq.per.sample.df.csv' ) ),
           quote = FALSE,
           row.names = TRUE )

# Print Cluster Frequencies in Stacked Bar Plot
file.name <- paste( section, 'freq', 'per', 'sample', 'pdf', sep='.' )

pdf( file = here( 'results', 'calculations', 'VDJ', file.name ) )

ggplot( data = counts, aes( fill=Cluster,
                            y=Counts,
                            x=factor(Sample, level = njoin) ) ) + # Samples are ordered by njoin
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual( values = sample(col_vector, n), aesthetics = "fill" ) + # Using 'n' defined above
  theme(panel.background=element_rect(fill="white"), # background = white
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"))

dev.off()

cat( paste0(section, " ", "section complete.\n") )
rm( section, n, counts, freq )


############################################################################################################
### Calculate DEGs for CSF-enriched clones compared to Blood-specific clones in the blood of NS subjects ###
############################################################################################################
# Load calculated CSF & PBMC frequencies / ratios for clonotypes
section <- 'Tissue.for.CD8.clones'
subject.names <- unique(MetaData$Subject)
nsubjects <- length( subject.names )

DefaultAssay(scd) <- "RNA"
scd <- NormalizeData(scd) # Normalize data
all.genes <- rownames(scd)
scd <- ScaleData(scd, features = all.genes) # Scale data

Idents(scd) <- 'lin_1'
tsub <- subset( scd, idents = 'CD8T' ) # Subset CD8 T cells

clone.freq <- list()
for ( i in 1:nsubjects ) {
  clone.file.name <- paste(
    'CD8' ,
    'clones' ,
    'CSF' ,
    'vs' ,
    'PBMC' ,
    subject.names[i] ,
    'csv' ,
    sep = '.' )

  clone.path <- here( 'results',
                      'calculations',
                      'VDJ',
                      clone.file.name )
  clone.freq[[i]] <- read.csv(clone.path)
  rm( clone.file.name, clone.path )
}

# Filter for CSF- and Blood-specific clonotypes
csf.clones <- list()
blood.clones <- list()
both.clones <- list()
for (i in 1:nsubjects ) {
  csf.clones[[i]] <- clone.freq[[i]] %>%
    filter( CSF_PBMC_ratio > 5 ) %>%
    filter( CSF > 5)

  blood.clones[[i]] <- clone.freq[[i]] %>%
    filter( CSF_PBMC_ratio < 0.2 ) %>%
    filter( PBMC > 5)

  both.clones[[i]] <- clone.freq[[i]] %>%
    filter( CSF_PBMC_ratio < 5 ) %>%
    filter( CSF_PBMC_ratio > 0.2 ) %>%
    filter( CSF > 0 ) %>%
    filter( PBMC > 0 )
}

# assign blood- vs CSF-enriched
tissue.specificity <- data.frame( barcodes = rownames( tsub@meta.data ),
                                  Subjects = tsub@meta.data$Subjects,
                                  tcr.seq = tsub@meta.data$tcr.seq )
rownames( tissue.specificity ) <- colnames( tsub )

tissue.specificity <- select( tissue.specificity, -barcodes )

tissue.specificity$clone.spec = 'neither'

for (i in 1:nrow(tissue.specificity) )
{
  for (j in 1:nsubjects)
  {
    if ( tissue.specificity$Subjects[i] == subject.names[j] & tissue.specificity$tcr.seq[i] %in% csf.clones[[j]]$tcr.seq ) {
      tissue.specificity$clone.spec[i] <- 'CSF-enriched'
    } else if ( tissue.specificity$Subjects[i] == subject.names[j] & tissue.specificity$tcr.seq[i] %in% blood.clones[[j]]$tcr.seq  ) {
      tissue.specificity$clone.spec[i] <- 'PBMC-specific'
    } else if ( tissue.specificity$Subjects[i] == subject.names[j] & tissue.specificity$tcr.seq[i] %in% both.clones[[j]]$tcr.seq  ) {
      tissue.specificity$clone.spec[i] <- 'non-specific'
    }
  }
}

table(tissue.specificity$clone.spec)

tsub <- AddMetaData( tsub,
                     tissue.specificity$clone.spec,
                     col.name = 'tissue.for.clone.2' )

##### ATTN here
# subset out blood CD8 T cells from NS subjects
Idents(tsub) <- 'Tissue'
blood.sub <- subset( tsub, idents = 'PBMC' )
Idents(blood.sub) <- 'Status'
blood.sub <- subset( blood.sub, idents = c( 'active', 'inactive' ) )

# subset out blood- vs CSF-enriched clonotypes
Idents(blood.sub) <- 'tissue.for.clone.2'
spec.sub <- subset( blood.sub, idents = c( 'CSF-enriched', 'PBMC-specific' ) )

# run FindAllMarkers and save DEGs
Idents(spec.sub) <- 'tissue.for.clone.2'
DefaultAssay(spec.sub) <- "RNA"
markers <- FindAllMarkers(spec.sub,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

x <- filter( markers, p_val_adj < 0.05 ) %>%
  group_by(cluster)

file.name <- paste( section, 'tsv', sep = '.' )
write.table( x = x,
             file = here( 'results', 'DEGs', file.name ),
             sep="\t",
             quote = FALSE,
             row.names = FALSE )

# print violin plots of markers for blood- and CSF-enriched clonotypes
Idents(tsub) <- 'Status'
clone.sub <- subset( tsub , idents = c( 'active', 'inactive' ) )
Idents(clone.sub) <- 'tissue.for.clone.2'
clone.sub <- subset( clone.sub, idents = c( 'CSF-enriched', 'PBMC-specific' ) )
DefaultAssay(clone.sub) <- 'RNA'

Idents(clone.sub) <- 'Tissue'
plot_list <- list()
for ( i in 1:nrow(x) ) {
  plot_list[[i]] <- VlnPlot(clone.sub,
                            features = x$gene[i],
                            cols = c('blue', 'red'),
                            pt.size = 0,
                            split.by = 'tissue.for.clone.2' )
}

for ( i in 1:nrow(x) ) {
  file.name <- paste( 'Violin',
                      'Only',
                      'Tissue',
                      'Specific',
                      'Clones',
                      x$gene[i],
                      'pdf',
                      sep = '.')

  pdf( file = here( 'results', 'DEGs', file.name ),
       width = 8,
       height = 11 )
  print( plot_list[[i]] )
  dev.off()
}

# Save the object
saveRDS(tsub,
        file = here( 'results' , 'objects' , 'CD8T.rds' ) )

############
# MONOCLE2 #
############

# Set named CD8T clusters (e.g. CD8T_1) as Idents
Idents(tsub) <- 'lin_2'
Idents(tsub) <- factor(x = Idents(tsub),
                       levels = sort(levels(tsub))) # Need to re-sort the Idents (which were factors that were out of order)
DefaultAssay(tsub) <- 'RNA'
# Convert Seurat Object to CellDataSet (cds) for Monocle2

# Extract data, phenotype data, and feature data from the SeuratObject
data <- GetAssayData(object = tsub,
                     slot = "counts") # REF: https://github.com/cole-trapnell-lab/monocle-release/issues/262
# data <- as(as.matrix(tsub@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = tsub@meta.data)

fData <- data.frame(gene_short_name = row.names(data),
                    row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

## Remove TCR (TRAV- & TRBV- from data (rows) and fd)
data <- data[ !grepl( "^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC|^IGHV|^IGHD|^IGHJ|^IGHM|^IGHE|^IGHG|^IGHA|^IGLV|^IGLD|^IGLC|^IGKV|^IGKD|^IGKC" , rownames(data) ) , ]
fd <- fd[ !grepl( "^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC|^IGHV|^IGHD|^IGHJ|^IGHM|^IGHE|^IGHG|^IGHA|^IGLV|^IGLD|^IGLC|^IGKV|^IGKD|^IGKC" , fd$gene_short_name ) , ]
data[1:5,1:5] # Sanity check
head(fd) # Sanity check

# intermediate cleanup
rm( list = c( 'tsub', 'scd', 'blood.sub' ) )
gc()

# Make CellDataSet object
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size() )
message("cds file created.")

# Estimate size factors and dispersions (REQUORED in Monocle2)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
message("Estimate Size Factorw and Dispersions complete.")

# Clustering cells without marker genes: http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
cds <- detectGenes(cds, min_expr = 0.1)
fData(cds)$use_for_ordering <- fData(cds)$num_cells_expressed > 0.05 * ncol(cds)

# Print Elbow plot of PCA
plot1 <- plot_pc_variance_explained(cds, return_all = FALSE)
file.name <- paste( 'ElbowPlot' , 'pdf' , sep = '.' )
pdf( file = here( 'results', 'Monocle2' , file.name ) )
print( plot1 )
dev.off()
message("Elbow Plot printed.")

# Project with tSNE to 2 dimensions
cds <- reduceDimension(cds,
                       max_components = 2,
                       norm_method = 'log',
                       num_dim = 4,
                       reduction_method = 'tSNE',
                       verbose = TRUE)
message("reduceDimension complete.")

# Cluster Cell
cds_expressed_genes <-  row.names( subset(fData(cds), num_cells_expressed >= 10) )
cds <- clusterCells(cds, num_clusters = 5, verbose = FALSE) # Select 5 cluster to correspond to Seurat Object
clustering_DEG_genes <- differentialGeneTest( cds[cds_expressed_genes,],
                                              fullModelFormulaStr = '~Cluster',
                                              cores = 1 )
message("Cluster Cells complete.")

# Order cells based on selected genes with high dispersion across cells
disp_table <- dispersionTable(cds)
cds_ordering_genes <- subset(disp_table,
                             mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit) %>%
  pull( gene_id )

cds <- setOrderingFilter(cds, ordering_genes = cds_ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
message("orderCells by dispersionTable complete.")

# Print Tranjectories
cd8.cluster.colors <- brewer.pal(10, 'Paired')
file.name <- paste('UMAP', 'lin_2','pdf', sep='.')
pdf( file = here( 'results', 'Monocle2' , 'Disp_Table', file.name ),
     width = 8,
     height = 11 )
plot_cell_trajectory(cds,
                     color_by = 'lin_2',
                     cell_size = 1,
                     show_branch_points = FALSE) +
  scale_color_manual( values = cd8.cluster.colors )
dev.off()

tissue.colors <- c('blue','red')
file.name <- paste('UMAP', 'Tissue','pdf', sep='.')
pdf( file = here( 'results', 'Monocle2' , 'Disp_Table', file.name ),
     width = 8,
     height = 11 )
plot_cell_trajectory(cds,
                     color_by = 'Tissue',
                     cell_size = 1,
                     show_branch_points = FALSE) +
  scale_color_manual( values = tissue.colors )
dev.off()

file.name <- paste('UMAP', 'Pseudotime','pdf', sep='.')
pdf( file = here( 'results', 'Monocle2' , 'Disp_Table', file.name ),
     width = 8,
     height = 11 )
plot_cell_trajectory(cds,
                     color_by = 'Pseudotime',
                     cell_size = 1,
                     show_branch_points = FALSE)
dev.off()

message("Trajectories by dispersionTable complete.")

# Save CDS file
saveRDS( cds, file = here('results', 'Monocle2' , 'CDS', 'CD8T.CDS.rds'))
message("New Monocle2 file saved.")

# Save individual data
save(data, pd, fData, fd, file = here('results', 'Monocle2' , 'CDS', 'CD8T.Rdata'))
message("Individual data object saved.")
