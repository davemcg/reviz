use std::io::{ BufRead, BufWriter, Write, stdout, stdin};
use std::fs::read_to_string;
use std::collections::{HashSet, HashMap};
use clap::{Arg, App};

static MSNG_CHAR:&str = "MSG";

fn main() {
    let matches = App::new("Parse Snaptron Counts")
                          .version("0.0.1")
                          .author("swamsaur")
                          .arg(Arg::with_name("whitelist")
                              .short("w")
                              .long("whitelist")
                              .required(false)
                              .takes_value(true)
                          )
                          .arg(Arg::with_name("genemapper")
                                .short("g")
                                .long("genemapper")
                                .required(false)
                                .takes_value(true)
                          )
                          .arg(Arg::with_name("sampleids")
                                .short("s")
                                .long("sampleids")
                                .required(false)
                                .takes_value(true)
                                )
                          .get_matches();
    
    // this is some real DOP product of italy level spaghetti 
    ////Input
    let mut whitelist_set:HashSet<String> = HashSet::with_capacity(800000);//potentially use all samples as a whitelist 
    let mut gene_mapper:HashMap<String, String> = HashMap::with_capacity(800000);
    let mut sampleid_mapper:HashMap<String, usize> = HashMap::with_capacity(800000);
    ////output
    let mut writer = BufWriter::new(stdout());
    let stdin = stdin();
    let reader = stdin.lock();


    if let Some(whitelist_file) = matches.value_of("whitelist"){
        let wl_samples = read_to_string(whitelist_file)
                                .expect("Bad Whitelist File");               
        for sample in wl_samples.trim().split("\n"){
            whitelist_set.insert(String::from(sample));
        }
    } 

    if let Some(gene_mapper_file) = matches.value_of("genemapper"){
        
        let lines = read_to_string(gene_mapper_file).expect("failed to read gene_mapper_file");
        
        for exon_id in lines.trim().split("\n"){
            //exon_id,gene_name
            let line_split = exon_id.split(",").collect::<Vec<&str>>();
            &gene_mapper.insert(String::from(line_split[0]), String::from(line_split[1]));

        }
    }
    
    if let Some(sampleid_mapper_file) = matches.value_of("sampleids"){
        for (i, sample_id) in read_to_string(sampleid_mapper_file).expect("failed to read sampleids_mapper_file").trim().split("\n").enumerate(){
            &sampleid_mapper.insert(String::from(sample_id), i);
        }
    }

    if (sampleid_mapper.len()==0 && gene_mapper.len() > 0) || (sampleid_mapper.len() > 0 && gene_mapper.len() == 0){
        panic!("Missing input for either sampleid_mapper or gene_mapper")
    }
    if gene_mapper.len() > 0{
        let sample_ids_all = &sampleid_mapper.keys().cloned().collect::<Vec<String>>().join(",");
        let header = format!("{},{},{}\n", "gene_name", "snaptron_id", sample_ids_all);
        &writer.write_all(header.as_bytes());
    }

    for line in reader.lines() {
        let mut clean_line:String;
        if gene_mapper.len() > 0{
            clean_line =  parse_line_genemapper(line.unwrap(),  &gene_mapper, &sampleid_mapper);
        } else {
            clean_line =  parse_line_whitelist(line.unwrap(),  &whitelist_set);
        }
        
        &writer.write_all(clean_line.as_bytes());
    }
    //&writer.write_all("out of loop".as_bytes());   
}

fn parse_line_whitelist(line: String, whitelist_set:&HashSet<String>) -> String{
    let line_data_raw:Vec<&str> = line.split('|').collect();
    let snaptron_id = String::from(line_data_raw[0]);
    let msg = String::from(MSNG_CHAR);
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
    let msg = String::from(MSNG_CHAR);
    empty_vec.retain(|x| x != &msg);
    let clean_string = empty_vec.join("\n");
    if clean_string == ""{ // line may not contain any whitelist samples
        return String::from("");
    } else {
        return clean_string + "\n";
    }

}

fn parse_line_genemapper(line: String, gene_mapper:&HashMap<String, String>, sampleid_mapper:&HashMap<String, usize>) -> String{
    let line_data_raw:Vec<&str> = line.split('|').collect();
    let snaptron_id = String::from(line_data_raw[0]);
    let try_gene_name = gene_mapper.get(&snaptron_id);
    let gene_name = match try_gene_name {
        Some(try_gene_name) => try_gene_name.clone(),
        None => snaptron_id.clone()
    };
    //line_data_raw[1] = sample_counts
    let line_vec = line_data_raw[1].trim_matches(',').split(',').collect::<Vec<&str>>();
    let mut empty_vec = vec!["0"; sampleid_mapper.len()];
    for count_raw in line_vec.iter(){
        let count_split = count_raw.split(':').collect::<Vec<&str>>();
        let sample_id = String::from(count_split[0]);
        if count_split.len() < 2{
            // this is en empty data point, ie sample=0
            continue;
        } else if count_split.len() > 2{
            panic!("a record had multiple values, we don't know how to handle this");
        } else{    
            let idx = sampleid_mapper.get(&sample_id).expect("sample id not in mapper");
            empty_vec[*idx] = count_split[1];
        }
    }
    
    
    let clean_string = empty_vec.join(",");
    if clean_string == ""{ // line may not contain any whitelist samples
        return String::from("");
    } else {
        return format!("{},{},{}\n", gene_name,snaptron_id,clean_string);
    }

}

