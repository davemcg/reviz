library(tidyverse)
eiad_sample_table <- read_tsv('/data/swamyvs/EiaD_build/sampleTable_2019-05-25_tissues.tab') %>% 
  rename(sample_accession = `#sample_accession`)
snaptron_sample_table <- data.table::fread('recount3_sra3h_samples.tsv') %>% as_tibble 
snaptron_gtex_sample_table <- data.table::fread('recount3_gtex_samples.tsv') %>% as_tibble
all_snaptron <- bind_rows(snaptron_sample_table %>% select(rail_id, run_acc=external_id),
                          snaptron_gtex_sample_table %>% select(rail_id, run_acc))
eiad_in_snaptron <- all_snaptron %>% filter(run_acc %in% eiad_sample_table$run_accession)

write(eiad_in_snaptron$rail_id, 'Eiad_in_snaptron_whitelist.txt', sep = '\n')
