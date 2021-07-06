suppressMessages(library(recount3))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(cowplot))
args = commandArgs(trailingOnly=TRUE)

# mess of a function which returns EITHER (type "gene" or "exon")
# 1. gene level normalized count tibble with user given meta added (along with recount3 meta)
# 2. for ONE gene, tibble of exon level info (same meta as 1 above)

r3_data_puller <- function(type = 'gene', 
                           study_accession,
                           organism = 'human',
                           meta,
                           gene
){
  
  
  # e.g. SRP163355
  # only one!
  if (length(study_accession) != 1){
    print('ONLY ONE STUDY YO')
    stop()
  }
  project_accession <- study_accession
  
  # metadata file
  # tsv where the first column is the run_accession and the second column is the metadata field 
  # e.g. genotype or age or tissue or whatever
  colnames(meta) <- c('name', 'Condition')
  
  # will use the first gene given if multiple for the exon return
  if (type == 'exon'){
    gene_exon <- gene[1]
  }
  
  # have to give species
  r3_projects <- available_projects(organism = organism)
  # grab project info 
  proj_info <- subset(
    r3_projects,
    project == project_accession & project_type == "data_sources"
  )
  
  #print('Pull exon level quant')
  if (type == 'exon'){
    rse_study_exon <- create_rse(proj_info, type = 'exon')
    
    ## exon
    assay(rse_study_exon, "counts") <- transform_counts(rse_study_exon)
    assays(rse_study_exon)$TPM <- recount::getTPM(rse_study_exon, length_var = "score")
    
    # pull count info
    gene_exon_tpm <- assays(rse_study_exon)$counts[
      rowRanges(rse_study_exon) %>% 
        as_tibble(rownames = 'id') %>% 
        filter(grepl(gene_exon, 
                     gene_name, 
                     ignore.case= TRUE), 
               transcript_type == 'protein_coding') %>% 
        pull(recount_exon_id) %>% 
        unique(),]
    # add metadata
    gene_exon_tpm_meta <- gene_exon_tpm %>% 
      as_tibble(rownames = 'recount_exon_id') %>% 
      left_join(rowRanges(rse_study_exon) %>% 
                  as_tibble(rownames = 'id'), 
                by = 'recount_exon_id') 
    
    row.names(gene_exon_tpm_meta) <- gsub('\\|', ' ', row.names(gene_exon_tpm_meta))
    return(gene_exon_tpm_meta)
  } else if (type == 'gene'){
    
    #print('Pull gene level quant')
    rse_study_gene <- create_rse(proj_info, type = 'gene')
    
    # create counts and TPM
    ## gene
    assay(rse_study_gene, "counts") <- transform_counts(rse_study_gene)
    assays(rse_study_gene)$TPM <- recount::getTPM(rse_study_gene, length_var = "score")
    
    gene_quant <- assays(rse_study_gene)$counts[rowRanges(rse_study_gene) %>% 
                                                  as_tibble(rownames = 'id') %>% 
                                                  filter(grepl(paste0(gene, collapse = '|'), 
                                                               gene_name, ignore.case= TRUE)) %>% 
                                                  pull(gene_id) %>% unique(), ]
    
    gene_quant_meta <- gene_quant %>%
      as_tibble(rownames = 'gene_id') %>% 
      pivot_longer(cols = contains('SRR')) %>% 
      left_join(rowRanges(rse_study_gene) %>% as_tibble(), by = 'gene_id') %>% 
      left_join(meta, by = 'name')
    
    return(gene_quant_meta)
  } else {
    print('Give either `gene` or `exon`')
    stop()
  }
}



