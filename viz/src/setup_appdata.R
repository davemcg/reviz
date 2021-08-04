library(shiny)
library(tidyverse)
library(ComplexHeatmap)
library(fst)
library(DBI)
library(RSQLite)
library(recount3)
library(glue)
library(parallel)
library(data.table)

### run these to make fst's for app, change file paths accordingly
args <- commandArgs(trailingOnly = T)
prefix <- args[1]
## Format for the new long sqlite file is `snaptron_id,sample_id,counts` (snaptron_id is a custom gene identifier)
## both snaptron_id and sample_id will need to be mapped to gene_name and SRA_accession accordingly
## theoretically this could be done upstream during the db (with little channge in run time) but right now this is where we're at

## make_EiaD_st 

EiaD_st <-  read_tsv('app_data/sampleTable_2019-05-25_tissues.tab') %>% 
  dplyr::rename(sample_accession = `#sample_accession`)
snaptron_st <- data.table::fread('app_data/recount3_sra3h_samples.tsv') %>% as_tibble 
snaptron_gtex_st <- data.table::fread('app_data/recount3_gtex_samples.tsv') %>% as_tibble
all_snaptron <- bind_rows(snaptron_st%>% select(rail_id, run_accession=external_id),
                          snaptron_gtex_st %>% select(rail_id, run_accession = run_acc))
write(all_snaptron$rail_id, file = 'app_data/recount3_sampleids.txt', sep = '\n')
eiad_in_snaptron <- EiaD_st %>% inner_join(all_snaptron) %>% mutate(rail_id = as.character(rail_id))
write_fst(eiad_in_snaptron, 'app_data/eiad_snaptron_sample_table.fst')
ALL_STUDY_ACCESSIONS <- unique(eiad_in_snaptron$study_accession)
ALL_SRS_ACCESSIONS <- unique(eiad_in_snaptron$sample_accession)
save(ALL_STUDY_ACCESSIONS, ALL_SRS_ACCESSIONS, file = 'app_data/samples_and_studies.Rdata')


## cache_Snaptron_proj_info 
r3_projects <- available_projects(organism = 'human' ) %>% as_tibble()
write_fst(r3_projects, 'app_data/recount3_proj_info.fst')


## make_snaptron_id_converters 
gene_mapper <- read_delim('app_data/snaptronid2gene_name.txt', delim = '|', col_names = c('snaptron_id', 'raw')) %>% 
  mutate(gene_id = str_split(raw, ':') %>% sapply(function(x) x[1]), 
         gene_name = str_split(raw, ':') %>% sapply(function(x) x[2]), 
         gene_type = str_split(raw, ':') %>% sapply(function(x) x[3]),
         snaptron_id = as.character(snaptron_id)) %>% 
  select(-raw)
## I'm not totally sure how to interpret the way they annotate exons, should revisit this at a later date.
exon_mapper <- read_delim('app_data/snaptronid2exon_name.txt', delim = '|', col_names = c('snaptron_id', 'raw')) %>% 
  mutate(left = str_split(raw, ':') %>% sapply(function(x) x[1]), 
         right = str_split(raw, ':') %>% sapply(function(x) x[2]),
         gene_id = left
  ) %>% 
  left_join(gene_mapper %>% select(-snaptron_id)) %>% 
  mutate(gene_name = replace(gene_name, is.na(gene_name), gene_id[is.na(gene_name)]),
         gene_type = replace(gene_name, is.na(gene_name), 'missing gene_name/unclear annotation') ) %>% 
  select(-raw, -left, -right)
gene_mapper <- gene_mapper %>% mutate(snaptron_id = as.character(snaptron_id))
exon_mapper <- exon_mapper %>% mutate(snaptron_id = as.character(snaptron_id))
write_fst(gene_mapper, 'app_data/snaptron_id_2_gene_name.fst')
write_fst(exon_mapper, 'app_data/snaptron_id_2_exon_name.fst')
write_csv(exon_mapper %>% select(snaptron_id, gene_name), 'app_data/snaptron_id_2_exon_name.csv', col_names = F)
ALL_GENE_NAMES <- unique(gene_mapper$gene_name)
ALL_EXON_NAMES <- unique(exon_mapper$gene_name)
save(ALL_GENE_NAMES, ALL_EXON_NAMES, file = 'app_data/valid_gene_exon_names.Rdata')



## make exon count fsts 

make_exon_count_fsts <- function(gene){
  
  t_exons <- filter(exon_mapper, gene_name == gene) %>% pull(snaptron_id) %>% paste(collapse = ", ")
  cmd_stem <- "sqlite3 {t_db_file} 'SELECT snaptron_id, samples FROM intron WHERE snaptron_id IN ({t_exons})' |../parse_sample_counts/target/release/parse_sample_counts "
  srav3_cmd = glue(cmd_stem, t_db_file = "sql_files/sra3vh_exons.sqlite", t_exons = t_exons)
  gtex_cmd <- glue(cmd_stem, t_db_file = "sql_files/gtex_exons.sqlite", t_exons = t_exons)
  
  srav3_dt <- fread(cmd = srav3_cmd,header = F ) %>% as_tibble 
  gtex_dt <- fread(cmd = gtex_cmd, header = F) %>% as_tibble
  full_data <- bind_rows(srav3_dt, gtex_dt)
  if(nrow(full_data) == 0) {
    write(gene, 'make_exon_fst2.log', append = T, sep='\n')
    return()
  }
  colnames(full_data) <- c('snaptron_id', 'sample_id', 'count' )
  full_data <- full_data %>% mutate(snaptron_id = as.character(snaptron_id), 
                                    sample_id = as.character(sample_id))
  outfile <-glue('app_data/exon_fst/{gene}_.fst')
  write_fst(full_data, outfile)
  return()
}

print('making FSTs')
mclapply(unique(exon_mapper$gene_name),function(g) try(make_exon_count_fsts(gene=g),silent = T), mc.cores = 22)


####set up example data
## make metadata
target_genes <- c('ABCA4', 'VSX1', 'CRX', 'NRL', 'RHO', 'GAPDH')
target_samples <- EiaD_st %>% filter(study_accession != 'SRP012682') %>%
  group_by('Tissue') %>% do(sample_n(., 100))
system('mkdir -p app_data/example_data/app_data')
eiad_in_snaptron_ss <-eiad_in_snaptron %>% 
  filter(sample_accession %in% target_samples$sample_accession )
write_fst(eiad_in_snaptron_ss, 'app_data/example_data/app_data/eiad_snaptron_sample_table.fst')
ALL_STUDY_ACCESSIONS <- unique(eiad_in_snaptron_ss$study_accession)
ALL_SRS_ACCESSIONS <- unique(eiad_in_snaptron_ss$sample_accession)
save(ALL_STUDY_ACCESSIONS, ALL_SRS_ACCESSIONS, file = 'app_data/example_data/app_data/samples_and_studies.Rdata')
gene_mapper_ss <- gene_mapper %>% filter(gene_name %in% target_genes)
exon_mapper_ss <- exon_mapper %>% filter(gene_name %in% target_genes)
write_fst(gene_mapper_ss, 'app_data/example_data/app_data/snaptron_id_2_gene_name.fst')
write_fst(exon_mapper_ss, 'app_data/example_data/app_data/snaptron_id_2_exon_name.fst')
write_csv(exon_mapper_ss %>% select(snaptron_id, gene_name), 'app_data/example_data/app_data/snaptron_id_2_exon_name.csv', col_names = F)
ALL_GENE_NAMES <- unique(gene_mapper_ss$gene_name)
ALL_EXON_NAMES <- unique(exon_mapper_ss$gene_name)
save(ALL_GENE_NAMES, ALL_EXON_NAMES, file = 'app_data/example_data/app_data/valid_gene_exon_names.Rdata')


## make subset sqlite data 

t_sample_id <- eiad_in_snaptron_ss %>% pull(rail_id)# haven't checked if snaptron_id<>gene_name is 1:1
t_snaptron_id <- gene_mapper_ss %>% pull(snaptron_id) %>% as.character()

prefix = 'recount3All'
db_file=glue("sql_files/{prefix}_gene_reformatted.db")
con<- dbConnect(SQLite(), dbname = db_file)
count_data <- tbl(con,'gene_counts_long') %>% 
  filter(snaptron_id %in% t_snaptron_id, sample_id %in% t_sample_id) %>%
  collect()
dbDisconnect(con)
system('mkdir -p app_data/example_data/sql_files')
ss_db_file <- glue("app_data/example_data/sql_files/{prefix}_gene_reformatted.db")
new_con <- dbConnect(RSQLite::SQLite(), ss_db_file)
dbWriteTable(new_con , 'gene_counts_long', count_data)
dbSendQuery(new_con , 'CREATE INDEX idx_gene ON gene_counts_long(snaptron_id, sample_id)')
dbDisconnect(new_con )

### make exon fsts 

system('mkdir -p app_data/example_data/app_data/exon_fst/')
make_exon_count_fsts_ss <- function(gene, t_sample_ids){
  
  t_exons <- filter(exon_mapper, gene_name == gene) %>% pull(snaptron_id) %>% paste(collapse = ", ")
  cmd_stem <- "sqlite3 {t_db_file} 'SELECT snaptron_id, samples FROM intron WHERE snaptron_id IN ({t_exons})' |../parse_sample_counts/target/release/parse_sample_counts "
  srav3_cmd = glue(cmd_stem, t_db_file = "sql_files/sra3vh_exons.sqlite", t_exons = t_exons)
  
  
  srav3_dt <- fread(cmd = srav3_cmd,header = F ) %>% as_tibble
  
  colnames(srav3_dt) <- c('snaptron_id', 'sample_id', 'counts' )
  outfile <-glue('app_data/example_data/app_data/exon_fst/{gene}_.fst')
  full_data <- srav3_dt %>% filter(sample_id %in% t_sample_ids) %>%
    mutate(snaptron_id = as.character(snaptron_id), 
           sample_id = as.character(sample_id))
    
  write_fst(full_data, outfile)
  return()
}
lapply(gene_mapper_ss$gene_name, function(x) make_exon_count_fsts_ss(x, t_sample_id))


