use std::io::{ BufRead, BufWriter, Write, stdout, stdin};
use std::fs::read_to_string;
use std::collections::HashSet;
use clap::{Arg, App};

fn main() {
    let matches = App::new("Parse Snaptron Counts")
                          .version("0.0.1")
                          .author("swamsaur")
                          .arg(Arg::with_name("whitelist")
                              .short("w")
                              .long("whitelist")
                              .required(false)
                              .takes_value(true)
                          ).get_matches();
    let mut whitelist_set:HashSet<String> = HashSet::with_capacity(800000);//potentially use all samples as a whitelist 
    if let Some(whitelist_file) = matches.value_of("whitelist"){
        let wl_samples = read_to_string(whitelist_file)
                                .expect("Bad Whitelist File");               
        for sample in wl_samples.split("\n"){
            whitelist_set.insert(String::from(sample));
        }
    } 

    let mut  writer = BufWriter::new(stdout());
    let stdin = stdin();
    let reader = stdin.lock();
    for line in reader.lines() {
        let clean_line =  parse_line(line.unwrap(),  &whitelist_set);
        &writer.write_all(clean_line.as_bytes());
    }
    //&writer.write_all("out of loop".as_bytes());   
}

fn parse_line(line: String, whitelist_set:&HashSet<String>) -> String{
    // raw line looks like: "<snaptron_id>|<sample_id>:<counts>"
    let line_data_raw:Vec<&str> = line.split('|').collect();
    let snaptron_id = String::from(line_data_raw[0]);
    let msg = String::from("N");
    //line_data_raw[1] = sample_counts
    let line_vec = line_data_raw[1].trim_matches(',').split(',').collect::<Vec<&str>>();
    let mut empty_vec = vec![msg; line_vec.len()];
    for (i, count_raw) in line_vec.iter().enumerate(){
        let count_split = count_raw.split(':').collect::<Vec<&str>>();
        let sample_id = String::from(count_split[0]);
        if count_split.len() < 2{
            // this is en empty data point, ie sample=0
            continue;
        } else if count_split.len() > 2{
            panic!("a record had multiple values, we don't know how to handle this");
        } else{    
            
            if whitelist_set.len() >0 && !whitelist_set.contains(&sample_id){
                continue
            }else{
            let formatted_counts = format!("{},{},{}", snaptron_id,count_split[0],count_split[1]);
            empty_vec[i] = formatted_counts; 
            }
        }
    }
    let msg = String::from("N");
    empty_vec.retain(|x| x != &msg);
    let clean_string = empty_vec.join("\n");
    return clean_string + "\n";
}