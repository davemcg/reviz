library(recount3)
library(tidyverse)
library(ComplexHeatmap)

# mouse lens (eye)
project_accession <- "SRP104440"
# violin plot genes
## pull gene data for a few lens related TF and crystallins (structural component of lens)
gene <- c('PAX6$', 'CRYGA', 'CRYB', 'CRYA', 'HSF4','PROX1$', 'RARB','SIX3$')

# heatmap gene 
gene_exon <- 'PAX6$'

r3_projects <- available_projects(organism = 'mouse')



proj_info <- subset(
  r3_projects,
  project == project_accession & project_type == "data_sources"
)

print('Pull exon level quant')
rse_study_exon <- create_rse(proj_info, type = 'exon')
# rse_study_exon

print('Pull gene level quant')
rse_study_gene <- create_rse(proj_info, type = 'gene')
# rse_study_gene

# create counts and TPM
## gene
assay(rse_study_gene, "counts") <- transform_counts(rse_study_gene)
assays(rse_study_gene)$TPM <- recount::getTPM(rse_study_gene, length_var = "score")
## exon
assay(rse_study_exon, "counts") <- transform_counts(rse_study_exon)
assays(rse_study_exon)$TPM <- recount::getTPM(rse_study_exon, length_var = "score")
#colSums(assay(rse_study_exon, "TPM")) / 1e6 ## Should all be equal to 1

# pull TPM quant
## broken for mouse for some reason
gene_quant <- assays(rse_study_gene)$counts[rowRanges(rse_study_gene) %>% 
                                              as_tibble(rownames = 'id') %>% 
                                              filter(grepl(paste0(gene, collapse = '|'), 
                                                           gene_name, ignore.case= TRUE)) %>% 
                                              pull(gene_id) %>% unique(), ]
###################################
# violin plot
###################################
meta <- colData(rse_study_gene) %>% 
  as_tibble() %>% 
  mutate(name = external_id,
         Genotype = str_extract(sra.sample_title, '\\(.*') %>% 
           gsub('\\(|\\)|\\d+','',.) %>% 
           trimws() %>% 
           tolower())

## make gene plot
gene_quant %>%
  as_tibble(rownames = 'gene_id') %>% 
  pivot_longer(cols = contains('SRR')) %>% 
  left_join(rowRanges(rse_study_gene) %>% as_tibble(), by = 'gene_id') %>% 
  left_join(meta, by = 'name') %>% 
  ggplot(aes(x=gene_name, y=log2(value+1), fill = Genotype)) +
  geom_violin(scale = 'width', position = position_dodge(width = 1)) +
  #scattermore::geom_scattermore(pointsize = 3, position = position_dodge(width = 1)) +
  #geom_boxplot(width = 0.3)+
  cowplot::theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values = pals::brewer.set1(n=12) %>% unname()) +
  xlab('Gene')

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
source('src/complex_heatmap_annotation_builders.R')
out_col = tx_col_df_maker(gene_exon_tpm_meta) 
out_row <- tx_row_df_maker(meta)

user_colors <- list(out_row$colors)
row_color_name <- 'Genotype'
names(user_colors) <- row_color_name
ha_column = HeatmapAnnotation(df = out_col$df, col = out_col$colors, 
                              show_legend = FALSE)
ha_row <- rowAnnotation(df = out_row$df,
                        col = user_colors,
                        name = "Genotype")
## make heatmap plot
Heatmap(t(log2(gene_exon_tpm+1)), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        name = 'log2(Gene Expression+1)',
        col = viridis::viridis(n=10), 
        top_annotation = ha_column, right_annotation = ha_row)


