use std::collections::HashMap;

use needletail::Sequence;

pub fn parse_reads(file_path: &str) -> HashMap<String, String> {
    let mut reads = needletail::parse_fastx_file(&file_path).expect("Invalid path to reads file");
    let mut sequences = HashMap::new();
    while let Some(record) = reads.next() {
        let seqrec = record.expect("Invalid read");
        let seq = String::from_utf8_lossy(seqrec.sequence()).to_string();
        // remove \n and spaces from seq
        let seq = seq.replace("\n", "").replace(" ", "");
        let id = String::from_utf8_lossy(seqrec.id()).to_string().split(" ").next().unwrap().to_string();
        sequences.insert(id, seq);
    }
    sequences
}

pub fn print_reads_keys(reads: &HashMap<String, String>) {
    for key in reads.keys() {
        println!("Read ID: {}", key);
    }
}