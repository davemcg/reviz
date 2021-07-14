#!/bin/bash 
module load SQLite3
mkdir -p sql_files
wget -O sql_files/genes.sqlite http://snaptron.cs.jhu.edu/data/srav3h/genes.sqlite
wget -O recount3_sra3h_samples.tsv http://snaptron.cs.jhu.edu/data/srav3h/samples.tsv
##make sure parse_sample_counts is built with `cargo build --release`
mkdir -p chunks
sqlite3 sql_files/genes.sqlite 'SELECT snaptron_id, samples FROM "intron" ' | parallel --blocksize 1500m --pipe 'parse_sample_counts/target/release/parse_sample_counts >> chunks/clean_data.tsv'
#### have to run this on sql command line 
# sqlite3 sql_files/genes_long.db
# create table intron_long(snaptron_id, sample_id, counts);
# .mode csv 
# .import chunks/test_cd.tsv intron_long;
# CREATE INDEX idx ON intron_long(snaptron_id, sample_id)