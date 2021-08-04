#!/bin/bash 

set -e 
## requires ~ 120GB + 24 cpus + 48 hours 
module load SQLite
mkdir -p sql_files
prefix=$1 # make sure this is a valid path
geneCounts="${prefix}_gene.csv"
db="sql_files/${prefix}_gene_reformatted.db"
##################
##Pull raw data ##
##################

# wget -O sql_files/sra3vh_genes.sqlite http://snaptron.cs.jhu.edu/data/srav3h/genes.sqlite #//81Gb 
# wget -O app_data/recount3_sra3h_samples.tsv http://snaptron.cs.jhu.edu/data/srav3h/samples.tsv 
# wget -O sql_files/sra3vh_exons.sqlite http://snaptron.cs.jhu.edu/data/srav3h/exons.sqlite #//879Gb

# wget -O sql_files/gtex_genes.sqlite http://snaptron.cs.jhu.edu/data/gtexv2/genes.sqlite
# wget -O sql_files/gtex_exons.sqlite http://snaptron.cs.jhu.edu/data/gtexv2/exons.sqlite
# wget -O app_data/recount3_gtex_samples.tsv http://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv

##########################################
##convert gene level data to long format## 
##########################################

rm -f $geneCounts
rm -f $db 

## make sure parse_sample_counts is built with `cargo build --release`
## SRA-gene
sqlite3 sql_files/sra3vh_genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 500m --pipe "../parse_sample_counts/target/release/parse_sample_counts  >> ${geneCounts}"
echo "Finished Reading sra3vh gene counts"

## Gtex-gene
sqlite3 sql_files/gtex_genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 500m --pipe "../parse_sample_counts/target/release/parse_sample_counts  >> ${geneCounts}"
echo "Finished Reading gtex gene counts"

## do not add semi colons after the `.` statements
echo "CREATE TABLE gene_counts_long(snaptron_id, sample_id, counts); 
.mode csv  
.import ${geneCounts} gene_counts_long
CREATE INDEX idx_gene ON gene_counts_long(snaptron_id, sample_id);" | sqlite3 -batch $db 
echo "Finished Writing gene db"

rm $geneCounts
#####################
##make app metadata## 
#####################
## since we are pulling just the gene ids, we only need to check the srav3h databse
sqlite3 sql_files/sra3vh_genes.sqlite "SELECT snaptron_id, right_annotated FROM intron" > app_data/snaptronid2gene_name.txt    
sqlite3 sql_files/sra3vh_exons.sqlite "SELECT snaptron_id, right_annotated FROM intron" > app_data/snaptronid2exon_name.txt

mkdir -p app_data/exon_fst/
module load R/4.0.3 
Rscript src/setup_appdata.R $prefix
