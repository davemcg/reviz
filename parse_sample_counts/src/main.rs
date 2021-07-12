use std::fs;
use std::io::{ BufRead, BufWriter, Write, stdout, stdin};
use std::collections::HashMap;
//use regex::Regex;
//use lazy_static::lazy_static;
use itertools::Itertools;
fn main() {
    //let sample_counts_file = fs::File::open(stdin()).expect("bad input");
    //let outfile = fs::File::create("parsed_counts.tsv").unwrap();
    let mut  writer = BufWriter::new(stdout());
    let snaptron_ids = fs::read_to_string("/data/swamyvs/reviz/ids.txt").unwrap();
    let snaptron_ids = snaptron_ids.split("\n").collect::<Vec<&str>>();
    let mut id_mapper:HashMap<&str, usize> = HashMap::new();
    for i in 0..snaptron_ids.len() {
        id_mapper.insert(&snaptron_ids[i], i);
    }
    let stdin = stdin();
    let reader = stdin.lock();
    for line in reader.lines() {
        let clean_line =  parse_line(line.unwrap(), &id_mapper);
        &writer.write_all(clean_line.as_bytes());
        &writer.write_all("\n".as_bytes());
    }   
}


fn parse_line(line: String, mapper: &HashMap<&str, usize>) -> String{
    // lazy_static! {
    //     static ref sample_pat: Regex = Regex::new(r"^\d+").unwrap();
    //     static ref count_pat: Regex = Regex::new(r"\d+$").unwrap();
    // }
    let line_vec =  line.trim_matches(',').split(',').collect::<Vec<&str>>();
    let sample_counts = line_vec.iter().map(|x| x.split(':').collect_tuple().unwrap()).collect::<Vec<(&str, &str)>>();
    let mut empty_vec = vec!["0"; mapper.len()];
    for sample in sample_counts{
        let idx = mapper.get(sample.0).expect("sample id missing from railid_list");
        empty_vec[*idx] = sample.1;
    }
    let clean_string = empty_vec.join(",");
    return clean_string;
}