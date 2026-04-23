use clap::Parser;
use rg_exp::{
    chain::{self, Chain},
    graph_index, parser, reads,
};
use std::{collections::HashMap, process::Command};
fn main() {
    let args = parser::Args::parse();
    let graph_path = args.graph;
    let reads_path = args.reads;

    let chain_bytes: Vec<u8> = if let Some(chain_path) = args.chain {
        eprintln!("Reading chain from file: {}", chain_path);
        std::fs::read(&chain_path).unwrap_or_else(|e| {
            eprintln!("Failed to read chain file: {}", e);
            std::process::exit(1);
        })
    } else {
        let output = Command::new("minigraph")
            .arg("-xlr")
            .arg("-c")
            .arg("-S")
            .arg("-t1")
            .arg(&graph_path)
            .arg(&reads_path)
            .output();
        match output {
            Ok(o) => {
                eprintln!("Minigraph command executed with status: {}", o.status);
                std::fs::write("minigraph_output.txt", &o.stdout)
                    .expect("Unable to write minigraph output to file");
                o.stdout
            }
            Err(e) => {
                eprintln!("Failed to execute minigraph: {}", e);
                std::process::exit(1);
            }
        }
    };

    // build path index
        let index = graph_index::index_paths(&graph_path);
        // load the reads
        let reads = reads::parse_reads(&reads_path);
        // split the output per read
        let split_alignments = split_reads_alignment(&chain_bytes);

        // gafs
        let mut gafs: HashMap<String, rg_exp::gaf::GAFStruct> = HashMap::new();

        // iterate over each read anchor
        for (read_id_b, chain_b) in split_alignments.iter() {
            // iterate over anchors, 1-2, 2-3, ...
            if chain_b.anchors.len() < 2 {
                continue;
            }
            eprintln!(
                "Processing read {} with {} anchors",
                String::from_utf8_lossy(&read_id_b),
                chain_b.anchors.len()
            );
            for window in chain_b.anchors.windows(2) {
                let from = &window[0];
                let to = &window[1];
                let mut rc = false;
                //println!("FROM: {:?}, TO: {:?}", from, to);
                // get read sequence slice delimited by from and to
                let mut read_slice: String;
                if let Some(read_seq) = reads.get(&String::from_utf8_lossy(&read_id_b).to_string())
                {
                    read_slice = read_seq[from.read_start..to.read_end].to_string();
                    if !from.graph_pos.orientation {
                        rc = true;
                        // reverse complement
                        let revcomp: String = read_slice
                            .chars()
                            .rev()
                            .map(|c| match c {
                                'A' => 'T',
                                'T' => 'A',
                                'C' => 'G',
                                'G' => 'C',
                                'N' => 'N',
                                _ => c,
                            })
                            .collect();
                        read_slice = revcomp;
                    }
                } else {
                    eprintln!(
                        "Read ID {} not found in reads",
                        String::from_utf8_lossy(&read_id_b)
                    );
                    read_slice = String::new();
                }
                // build the subgraph between from and to
                let subgraph = chain::extract_subgraph(from, to, &index, read_slice.len());
                let mut gfa_out = String::new();
                gfa::writer::write_gfa(&subgraph, &mut gfa_out);
                gfa_out.insert_str(0, "GRAPH:\n");
                let read_out = format!(
                    "READ:\n>{}\n{}",
                    String::from_utf8_lossy(&read_id_b),
                    read_slice
                );
                // DEBUG
                eprintln!("Subgraph GFA:\n{}", gfa_out);
                eprintln!("Read slice:\n{}", read_out);

                // do alignment with recalign
                let run_recalign = Command::new("../recalign/target/release/recalign")
                    .arg("-efast")
                    .arg("-s8")
                    .arg("-m")
                    .arg("-k2")
                    .arg("-r1")
                    .arg("-g")
                    .arg("-") // read graph from stdin
                    .arg("-q")
                    .arg("-") // read query from stdin
                    .stdin(std::process::Stdio::piped())
                    .stdout(std::process::Stdio::piped())
                    .spawn()
                    .and_then(|mut child| {
                        {
                            let stdin = child.stdin.as_mut().expect("Failed to open stdin");
                            use std::io::Write;
                            stdin
                                .write_all(gfa_out.as_bytes())
                                .expect("Failed to write graph to stdin");

                            stdin
                                .write_all(read_out.as_bytes())
                                .expect("Failed to write read slice to stdin");
                        }
                        child.wait_with_output()
                    });

                match run_recalign {
                    Ok(recalign_output) => {
                        eprintln!("Recalign executed with status: {}", recalign_output.status);
                        if !recalign_output.status.success() {
                            eprintln!(
                                "Recalign failed for anchor pair ({} -> {}), skipping.",
                                String::from_utf8_lossy(&from.graph_pos.node_id),
                                String::from_utf8_lossy(&to.graph_pos.node_id)
                            );
                            continue;
                        }

                        let mut gaf = get_gaf_from_recalign_output(&recalign_output.stdout, rc);
                        // fix: use actual read length, not the slice length recalign saw
                        gaf.query_length = reads
                            .get(&String::from_utf8_lossy(&read_id_b).to_string())
                            .map(|s| s.len())
                            .unwrap_or(0);
                        // update read_start and read_end to be relative to the whole read
                        gaf.query_start += from.read_start;
                        gaf.query_end += from.read_start;

                        // if first time, insert, else, merge
                        let read_id_str = String::from_utf8_lossy(&read_id_b).to_string();
                        if let Some(existing_gaf) = gafs.get_mut(&read_id_str) {
                            let merged_gaf =
                                rg_exp::gaf::GAFStruct::merge_gafs(existing_gaf, &gaf);
                            *existing_gaf = merged_gaf;
                        } else {
                            gafs.insert(read_id_str, gaf);
                        }
                    }
                    Err(e) => {
                        eprintln!("Failed to execute recalign: {}", e);
                    }
                }
            }
        }
        // print gafs
        for (read_id, gaf) in gafs.iter() {
            println!("GAF for read {}:\n{}", read_id, gaf.clone().to_string());
        }
}

fn split_reads_alignment(output: &[u8]) -> Vec<(Vec<u8>, Chain)> {
    let mut reads_alignments = Vec::new();
    let mut current_read_alignment = String::new();
    let mut read_id: Vec<u8> = Vec::new();
    let lines = output.split(|&b| b == b'\n').into_iter();
    //first line is always the read ID
    for line in lines {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'*' {
            current_read_alignment.push_str(&String::from_utf8_lossy(line));
            current_read_alignment.push('\n');
        } else {
            if !current_read_alignment.is_empty() {
                let chain = chain::build_chain(&current_read_alignment);
                reads_alignments.push((read_id.clone(), chain));
                current_read_alignment.clear();
            }
            read_id = line.split(|&b| b == b'\t').next().unwrap().to_vec();
        }
    }
    // flush last read
    if !current_read_alignment.is_empty() {
        let chain = chain::build_chain(&current_read_alignment);
        reads_alignments.push((read_id.clone(), chain));
    }
    /*
    println!(
        "Total reads alignments parsed: {:#?}",
        reads_alignments.len()
    );
    for (rid, chain) in reads_alignments.iter() {
        println!(
            "Read ID: {}, Alignment:\n{:#?}",
            String::from_utf8_lossy(&rid),
            chain
        );
    }
     */
    reads_alignments
}

fn get_gaf_from_recalign_output(output: &[u8], rc: bool) -> rg_exp::gaf::GAFStruct {
    let mut lines = output.split(|&b| b == b'\n').into_iter();
    // two lines, concatenate them and parse
    let mut gaf_line = String::new();
    gaf_line.push_str(&String::from_utf8_lossy(lines.next().unwrap()));
    // for second line, it's only comments, change \t to , with extra \t at start 
    let comments_line = lines.next().unwrap();
    if !comments_line.is_empty() {
        gaf_line.push('\t');
        let comments_str = String::from_utf8_lossy(comments_line).replace("\t", ",");
        gaf_line.push_str(&comments_str);
    }
    let mut gaf = rg_exp::gaf::GAFStruct::from_string(&gaf_line);
    if rc {
        gaf = gaf.reverse_complement();
    }
    gaf
}
