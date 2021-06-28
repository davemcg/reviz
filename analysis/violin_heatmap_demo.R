library(recount3)
library(tidyverse)
library(ComplexHeatmap)
r3_projects <- available_projects(organism = 'mouse')

proj_info <- subset(
  r3_projects,
  project == "SRP104440" & project_type == "data_sources"
)


rse_study_exon <- create_rse(proj_info, type = 'exon')
rse_study_exon

rse_study_gene <- create_rse(proj_info, type = 'gene')
rse_study_gene
# just raw counts
# assays(rse_study_exon)
# id ranges for a gene of interest
# rowRanges(rse_study_exon) %>% as_tibble() %>% filter(grepl('PAX6', gene_name, ignore.case= TRUE)) %>% group_by(transcript_id) %>% count()

# create counts and TPM
## gene
assay(rse_study_gene, "counts") <- transform_counts(rse_study_gene)
assays(rse_study_gene)$TPM <- recount::getTPM(rse_study_gene, length_var = "score")
## exon
assay(rse_study_exon, "counts") <- transform_counts(rse_study_exon)
assays(rse_study_exon)$TPM <- recount::getTPM(rse_study_exon, length_var = "score")
#colSums(assay(rse_study_exon, "TPM")) / 1e6 ## Should all be equal to 1

# pull gene data for just pax6 and abca4
gene <- c('PAX6$', 'CRYGA', 'CRYB', 'CRYA', 'HSF4','PROX1$', 'RARB','SIX3$')
gene_tpm <- assays(rse_study_gene)$TPM[rowRanges(rse_study_gene) %>% 
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

gene_tpm %>%
  as_tibble(rownames = 'gene_id') %>% 
  pivot_longer(cols = contains('SRR')) %>% 
  left_join(rowRanges(rse_study_gene) %>% as_tibble(), by = 'gene_id') %>% 
  left_join(meta, by = 'name') %>% 
  ggplot(aes(x=gene_name, y=log2(value+1), color = Genotype)) +
  geom_violin(scale = 'width') +
  #geom_boxplot(width = 0.3)+
  cowplot::theme_cowplot()

#################
# exon heatmap
#################

# pull exon data for just pax6 exons
gene <- 'CRYAA$'
gene_exon_tpm <- assays(rse_study_exon)$counts[rowRanges(rse_study_exon) %>% as_tibble(rownames = 'id') %>% filter(grepl(gene, gene_name, ignore.case= TRUE), transcript_type == 'protein_coding') %>% data.frame() %>% pull(recount_exon_id) %>% unique(),]
# add metadata
gene_exon_tpm_meta <- gene_exon_tpm %>% as_tibble(rownames = 'recount_exon_id') %>% left_join(rowRanges(rse_study_exon) %>% as_tibble(rownames = 'id'), by = 'recount_exon_id')
# heatmap
## col metadata annotation
### arbitrary tx anno df function
tx_col_df_maker <- function(meta){
  df <- gene_exon_tpm_meta %>% 
    select(transcript_id, recount_exon_id) %>% 
    unique()
  for (i in unique(df$transcript_id)){
    df[,i] <- grepl(i, df$transcript_id)
  }
  out <- list()
  out$df <- df %>% select(-transcript_id, -recount_exon_id) %>% data.frame()
  out$df[out$df == TRUE] <- 'TRUE'
  out$df[out$df == FALSE] <- 'FALSE'
  colors <- list()
  for (i in colnames(df)){
    colors[[i]] <- c('TRUE' = 'black', 'FALSE' = 'white')
  }
  out$colors <- colors
  out
}

tx_row_df_maker <- function(meta){
  df <- data.frame(Genotype = t(gene_exon_tpm) %>% as_tibble(rownames = 'name') %>% 
                     left_join(meta, by = 'name') %>% pull(Genotype))
  colors <- c(pals::alphabet(n=26), pals::alphabet2(n=26), pals::glasbey(n=32))
  colors <- rep(colors, 100)
  names(colors) <- unique(df[,1])
  colors <- colors[!is.na(names(colors))]
  out <- list()
  out[['df']] <- df
  out[['colors']] <- colors
  out
  
  
}

out_col = tx_col_df_maker(gene_exon_tpm_meta) 
out_row <- tx_row_df_maker(meta)

user_colors <- list(out_row$colors)
row_color_name <- 'Genotype'
names(user_colors) <- row_color_name
ha_column = HeatmapAnnotation(df = out_col$df, col = out_col$colors, 
                              show_legend = FALSE, 
                              name = 'log2(TPM+1)' )
ha_row <- rowAnnotation(df = out_row$df,
                        col = user_colors,
                        name = "Genotype")
Heatmap(t(log2(gene_exon_tpm+1)), 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = viridis::viridis(n=10), 
        top_annotation = ha_column, right_annotation = ha_row)


