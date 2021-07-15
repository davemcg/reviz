#!/bin/bash 

## requires ~ 20GB + 8 cpus
module load SQLite
mkdir -p sql_files
prefix=$1 # make sure this is a valid path
geneCounts="${prefix}_gene.csv"
exonCounts="${prefix}_exon.csv"
db="${prefix}_reformatted.db"
# wget -O sql_files/sra3vh_genes.sqlite http://snaptron.cs.jhu.edu/data/srav3h/genes.sqlite //81Gb 
# wget -O recount3_sra3h_samples.tsv http://snaptron.cs.jhu.edu/data/srav3h/samples.tsv
# wget -O sql_files/sra3vh_exons.sqlite http://snaptron.cs.jhu.edu/data/srav3h/exons.sqlite //879Gb

# wget -O sql_files/gtex_genes.sqlite http://snaptron.cs.jhu.edu/data/gtexv2/genes.sqlite
# wget -O sql_files/gtex_exons.sqlite http://snaptron.cs.jhu.edu/data/gtexv2/exons.sqlite
# wget -O recount3_gtex_samples.tsv http://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv

## make sure parse_sample_counts is built with `cargo build --release`
## run src/make_snaptron_whitelist.R to make whietlist
rm -f $geneCounts
rm -f $exonCounts
rm -f $db 
## SRA-gene
sqlite3 sql_files/sra3vh_genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 500m --pipe "parse_sample_counts/target/release/parse_sample_counts -w Eiad_in_snaptron_whitelist.txt >> ${geneCounts}"
echo "Finished Reading sra3vh gene counts"

## SRA-exon
sqlite3 sql_files/sra3vh_exons.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 500m --pipe "parse_sample_counts/target/release/parse_sample_counts -w Eiad_in_snaptron_whitelist.txt >> ${exonCounts}"
echo "Finished Reading sra3vh exon counts"

## Gtex-gene
sqlite3 sql_files/gtex_genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 500m --pipe "parse_sample_counts/target/release/parse_sample_counts -w Eiad_in_snaptron_whitelist.txt >> ${geneCounts}"
echo "Finished Reading gtex gene counts"

## Gtex-exon
sqlite3 sql_files/gtex_exons.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 500m --pipe "parse_sample_counts/target/release/parse_sample_counts -w Eiad_in_snaptron_whitelist.txt >> ${exonCounts}"
echo "Finished Reading gtex exon counts"

echo "CREATE TABLE gene_counts_long(snaptron_id, sample_id, counts); 
.mode csv ; 
.import ${geneCounts} gene_counts_long; 
CREATE INDEX idx_gene ON gene_counts_long(snaptron_id, sample_id);" | sqlite3 -batch $db 
echo "Finished Writing gene db"

echo "CREATE TABLE exon_counts_long(snaptron_id, sample_id, counts); 
.mode csv ; 
.import ${exonCounts} exon_counts_long; 
CREATE INDEX idx_exon ON exon_counts_long(snaptron_id, sample_id);" | sqlite3 -batch $db 
echo "Finished Writing exon db"