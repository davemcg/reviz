#!/bin/bash 

## requires ~ 20GB + 8 cpus
module load SQLite
mkdir -p sql_files
counts=$1
db=$2
# wget -O sql_files/sra3vh_genes.sqlite http://snaptron.cs.jhu.edu/data/srav3h/genes.sqlite //81Gb 
# wget -O recount3_sra3h_samples.tsv http://snaptron.cs.jhu.edu/data/srav3h/samples.tsv
# wget -O sql_files/sra3vh_exons.sqlite http://snaptron.cs.jhu.edu/data/srav3h/exons.sqlite //879Gb

# wget -O sql_files/gtex_genes.sqlite http://snaptron.cs.jhu.edu/data/gtexv2/genes.sqlite
# wget -O sql_files/gtex_exons.sqlite http://snaptron.cs.jhu.edu/data/gtexv2/exons.sqlite
# wget -O recount3_gtex_samples.tsv http://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv

## make sure parse_sample_counts is built with `cargo build --release`
## run src/make_snaptron_whitelist.R to make whietlist
rm -f $counts
rm -f $db 
mkdir -p counts_long

sqlite3 sql_files/sra3vh_genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 1000m --pipe "parse_sample_counts/target/release/parse_sample_counts -w Eiad_in_snaptron_whitelist.txt >> ${counts}"
echo "Finished Reading sra3vh"

sqlite3 sql_files/gtex_genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 1000m --pipe "parse_sample_counts/target/release/parse_sample_counts -w Eiad_in_snaptron_whitelist.txt >> ${counts}"

echo "Finished Reading GTEX"

echo "CREATE TABLE counts_long(snaptron_id, sample_id, counts); 
.mode csv ; 
.import ${counts} counts_long; 
CREATE INDEX idx ON counts_long(snaptron_id, sample_id);" | sqlite3 -batch $db 
echo "Finished Writing db"