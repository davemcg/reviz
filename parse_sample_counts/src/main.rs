use std::fs;
use std::io::{ BufRead, BufWriter, Write, stdout, stdin};
use std::collections::HashMap;
//use regex::Regex;
//use lazy_static::lazy_static;
//use itertools::Itertools;
fn main() {
    //let sample_counts_file = fs::File::open(stdin()).expect("bad input");
    //let outfile = fs::File::create("parsed_counts.tsv").unwrap();
    let mut  writer = BufWriter::new(stdout());
    let sample_ids = fs::read_to_string("/data/swamyvs/reviz/ids.txt").unwrap();
    let sample_ids = sample_ids.split("\n").collect::<Vec<&str>>();
    let mut id_mapper:HashMap<&str, usize> = HashMap::new();
    for i in 0..sample_ids.len() {
        id_mapper.insert(&sample_ids[i], i);
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
    let line_data_raw:Vec<&str> = line.split('|').collect();
    let snaptron_id = line_data_raw[0];
    let msg = String::from("N");
    //line_data_raw[1] = sample_counts
    let line_vec = line_data_raw[1].trim_matches(',').split(',').collect::<Vec<&str>>();
    let mut empty_vec = vec![msg; mapper.len()];
    for count_raw in line_vec{
        let count_split = count_raw.split(':').collect::<Vec<&str>>();
        if count_split.len() < 2{
            continue;
        } else if count_split.len() > 2{
            panic!("a record had multiple values, we don't know how to handle this");
        } else{    
            let idx = mapper.get(count_split[0]).expect("Sample id missing from rail_id metadata");
            let formatted_counts = format!("{},{}", snaptron_id, count_split[1]);
            empty_vec[*idx] = formatted_counts; 
        }
    }
    let msg = String::from("N");
    empty_vec.retain(|x| x != &msg);
    let clean_string = empty_vec.join("\n");
    return clean_string;
}