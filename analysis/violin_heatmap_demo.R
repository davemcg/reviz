suppressMessages(library(recount3))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(cowplot))
args = commandArgs(trailingOnly=TRUE)

# e.g. SRP163355
# only one!
project_accession <- args[1] #c("SRP163355")#,"SRP119291", "SRP091605")

# metadata file
# tsv where the first column is the run_accession and the second column is the metadata field 
# e.g. genotype or age or tissue or whatever
# NO COLUMN NAMES
meta <- read_tsv(args[2])
colnames(meta) <- c('name', 'Condition')
# violin/boxplot plot genes
## below genes are for a few lens related TF and crystallins (structural component of lens)
## gene <- c('PAX6$', 'CRYGA', 'CRYB', 'CRYA', 'HSF4','PROX1$', 'RARB','SIX3$', 'PRX$')
## give as one gene per line. regex's allowed.
## FIRST GENE will be the "exon" gene
gene <- scan(args[3], what = 'character')
# heatmap gene 
#gene_exon <- '^PRX$'
gene_exon <- gene[1]

# have to give species
r3_projects <- available_projects(organism = args[4])
# grab project info 
proj_info <- subset(
  r3_projects,
  project == project_accession & project_type == "data_sources"
)

#print('Pull exon level quant')
rse_study_exon <- create_rse(proj_info, type = 'exon')


#print('Pull gene level quant')
rse_study_gene <- create_rse(proj_info, type = 'gene')


# create counts and TPM
## gene
assay(rse_study_gene, "counts") <- transform_counts(rse_study_gene)
assays(rse_study_gene)$TPM <- recount::getTPM(rse_study_gene, length_var = "score")
## exon
assay(rse_study_exon, "counts") <- transform_counts(rse_study_exon)
assays(rse_study_exon)$TPM <- recount::getTPM(rse_study_exon, length_var = "score")
#colSums(assay(rse_study_exon, "TPM")) / 1e6 ## Should all be equal to 1

# pull TPM quant
####### broken for mouse for some reason #######
####### so just grabbing counts for now ########
gene_quant <- assays(rse_study_gene)$counts[rowRanges(rse_study_gene) %>% 
                                              as_tibble(rownames = 'id') %>% 
                                              filter(grepl(paste0(gene, collapse = '|'), 
                                                           gene_name, ignore.case= TRUE)) %>% 
                                              pull(gene_id) %>% unique(), ]
###################################
# violin plot
###################################

## custom bit
## to extract relevant sample info
# meta <- colData(rse_study_gene) %>% 
#   as_tibble() %>% 
#   mutate(name = external_id,
#          Condition = gsub('-\\d+','',sra.experiment_title)) %>% 
#   select(name, Condition)


## make gene plot
pdf(width = 8, height = 10, file = paste0(args[5], '_gene.pdf'))
gene_quant %>%
  as_tibble(rownames = 'gene_id') %>% 
  pivot_longer(cols = contains('SRR')) %>% 
  left_join(rowRanges(rse_study_gene) %>% as_tibble(), by = 'gene_id') %>% 
  left_join(meta, by = 'name') %>% 
  ggplot(aes(x=Condition, y=log2(value+1), fill = Condition)) +
  #geom_point() + 
  geom_boxplot() +
  #scattermore::geom_scattermore(pointsize = 3, position = position_dodge(width = 1)) +
  #geom_boxplot(width = 0.3)+
  theme_cowplot() +
  guides(x =  guide_axis(angle = 90)) +
  scale_fill_manual(values = c(pals::brewer.set1(n=12), pals::brewer.set2(n=10)) %>% unname()) +
  xlab('Gene') + facet_wrap(~gene_name, ncol = 5, scales = 'free_x')
dev.off()
#################
# exon heatmap
#################


# filter to only protein coding tx
gene_exon_tpm <- assays(rse_study_exon)$counts[
  rowRanges(rse_study_exon) %>% 
    as_tibble(rownames = 'id') %>% 
    filter(grepl(gene_exon, 
                 gene_name, 
                 ignore.case= TRUE), 
           transcript_type == 'protein_coding') %>% 
    pull(recount_exon_id) %>% 
    unique(),]
# grab metadata
gene_exon_tpm_meta <- gene_exon_tpm %>% 
  as_tibble(rownames = 'recount_exon_id') %>% 
  left_join(rowRanges(rse_study_exon) %>% 
              as_tibble(rownames = 'id'), 
            by = 'recount_exon_id') 

# build tx column info
source('../src/complex_heatmap_annotation_builders.R')
out_col = tx_col_df_maker(gene_exon_tpm_meta) 
out_row <- tx_row_df_maker(meta)

user_colors <- list(out_row$colors)
row_color_name <- 'Condition'
names(user_colors) <- row_color_name
ha_column = HeatmapAnnotation(df = out_col$df, col = out_col$colors, 
                              show_legend = FALSE)
ha_row <- rowAnnotation(df = out_row$df,
                        col = user_colors,
                        name = "Condition")
## make heatmap plot
pdf(width = 8, height = 10, file = paste0(args[5], '_exon.pdf'))
row.names(gene_exon_tpm) <- gsub('\\|', ' ', row.names(gene_exon_tpm))
Heatmap(t(log2(gene_exon_tpm+1)), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        name = 'log2(Gene Expression+1)',
        col = viridis::viridis(n=10), 
        top_annotation = ha_column,
        right_annotation = ha_row)
dev.off()


