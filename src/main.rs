use clap::Parser;
use rayon::prelude::*;
use rg_exp::{
    chain::{self, Chain},
    graph_index, parser, reads,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process::Command;

const RECALIGN_BIN: &str = "recalign/target/release/recalign";
const ALIGNMENT_CHUNK_SIZE: usize = 100;

fn main() {
    let args = parser::Args::parse();
    let graph_path = args.graph;
    let reads_path = args.reads;
    let max_reads = args.max_reads.max(1);

    // build path index
    let index = graph_index::index_paths(&graph_path);
    // load the reads
    let reads = reads::parse_reads(&reads_path);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(max_reads)
        .build()
        .expect("Failed to build rayon thread pool");

    if let Some(chain_path) = args.chain {
        eprintln!("Reading chain from file: {}", chain_path);
        let file = File::open(&chain_path).unwrap_or_else(|e| {
            eprintln!("Failed to open chain file {}: {}", chain_path, e);
            std::process::exit(1);
        });
        let reader = BufReader::new(file);
        process_chain_stream(reader, &reads, &index, &pool, None);
    } else {
        let mut child = Command::new("minigraph")
            .arg("-xlr")
            .arg("--vc")
            .arg("-c")
            .arg("-S")
            .arg("-N0")
            .arg("-t1")
            .arg(&graph_path)
            .arg(&reads_path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .unwrap_or_else(|e| {
                eprintln!("Failed to execute minigraph: {}", e);
                std::process::exit(1);
            });

        let stdout = child.stdout.take().unwrap_or_else(|| {
            eprintln!("Failed to capture minigraph stdout");
            std::process::exit(1);
        });
        let reader = BufReader::new(stdout);
        let mut out_file = File::create("minigraph_output.txt").unwrap_or_else(|e| {
            eprintln!("Unable to create minigraph_output.txt: {}", e);
            std::process::exit(1);
        });

        process_chain_stream(reader, &reads, &index, &pool, Some(&mut out_file));

        let status = child.wait().unwrap_or_else(|e| {
            eprintln!("Failed waiting for minigraph: {}", e);
            std::process::exit(1);
        });
        eprintln!("Minigraph command executed with status: {}", status);
    }
}

fn process_chain_stream<R: BufRead>(
    reader: R,
    reads: &HashMap<String, String>,
    index: &graph_index::PathIndex,
    pool: &rayon::ThreadPool,
    mut mirror_out: Option<&mut File>,
) {
    let mut pending: Vec<(Vec<u8>, Chain, String)> = Vec::with_capacity(ALIGNMENT_CHUNK_SIZE);
    let mut current_read_alignment = String::new();
    let mut current_primary_line = String::new();
    let mut read_id: Vec<u8> = Vec::new();
    let mut chunk_no: usize = 0;

    for line_res in reader.lines() {
        let line = line_res.unwrap_or_else(|e| {
            eprintln!("Failed to read alignment line: {}", e);
            std::process::exit(1);
        });

        if let Some(f) = mirror_out.as_mut() {
            if let Err(e) = writeln!(*f, "{}", line) {
                eprintln!("Failed writing minigraph output mirror: {}", e);
                std::process::exit(1);
            }
        }

        if line.is_empty() {
            continue;
        }
        if line.as_bytes()[0] == b'*' {
            current_read_alignment.push_str(&line);
            current_read_alignment.push('\n');
        } else {
            if !current_read_alignment.is_empty() {
                let chain = chain::build_chain(&current_read_alignment);
                pending.push((read_id.clone(), chain, current_primary_line.clone()));
                current_read_alignment.clear();
            }
            read_id = line
                .split('\t')
                .next()
                .unwrap_or_default()
                .as_bytes()
                .to_vec();
            current_primary_line = line;
        }

        if pending.len() >= ALIGNMENT_CHUNK_SIZE {
            chunk_no += 1;
            process_one_chunk(&pending, reads, index, pool, chunk_no);
            // Ensure any mirrored minigraph output is flushed to disk promptly.
            if let Some(f) = mirror_out.as_mut() {
                if let Err(e) = f.flush() {
                    eprintln!("Warning: failed to flush mirror file after chunk {}: {}", chunk_no, e);
                }
            }
            pending.clear();
        }
    }

    if !current_read_alignment.is_empty() {
        let chain = chain::build_chain(&current_read_alignment);
        pending.push((read_id.clone(), chain, current_primary_line.clone()));
    }

    if !pending.is_empty() {
        chunk_no += 1;
        process_one_chunk(&pending, reads, index, pool, chunk_no);
        if let Some(f) = mirror_out.as_mut() {
            if let Err(e) = f.flush() {
                eprintln!("Warning: failed to flush mirror file after final chunk {}: {}", chunk_no, e);
            }
        }
    }
}

fn process_one_chunk(
    chunk: &[(Vec<u8>, Chain, String)],
    reads: &HashMap<String, String>,
    index: &graph_index::PathIndex,
    pool: &rayon::ThreadPool,
    chunk_no: usize,
) {
    let mut gafs = pool.install(|| {
        chunk
            .par_iter()
            .enumerate()
            .filter_map(|(i, (read_id_b, chain_b, primary_line))| {
                process_chain_alignment(i, read_id_b, chain_b, primary_line, reads, index)
            })
            .collect::<Vec<_>>()
    });

    // Preserve deterministic output ordering inside each chunk.
    gafs.sort_by_key(|(i, _, _)| *i);

    // print gafs (one line per chain; multi-mapping reads yield multiple lines)
    for (_, _read_id, gaf) in gafs.iter() {
        println!("{}", gaf.clone().to_string());
    }

    // Force progressive writes even when stdout is redirected to file/pipe.
    // Lock stdout before flushing for more reliable behavior under heavy concurrency.
    {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        if let Err(e) = handle.flush() {
            eprintln!("Warning: failed to flush stdout after chunk {}: {}", chunk_no, e);
        }
    }
    eprintln!("Chunk {} flushed ({} chains).", chunk_no, chunk.len());
}

fn make_anchor_gaf(
    anchor: &rg_exp::chain::Anchor,
    read_id: &str,
    full_read_len: usize,
    path_offset: usize,
) -> rg_exp::gaf::GAFStruct {
    let node_id_num = std::str::from_utf8(&anchor.graph_pos.node_id)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let q_span = anchor.read_end.saturating_sub(anchor.read_start);
    let r_span = anchor.graph_pos.position.1.saturating_sub(anchor.graph_pos.position.0);
    let match_len = q_span.min(r_span);
    let comments = format!("tp:A:P\tNM:i:0\tdv:f:0.0000\tcg:Z:{}=", match_len);
    let path_dir = if anchor.graph_pos.orientation { '>' } else { '<' };
    rg_exp::gaf::GAFStruct::build_gaf_struct(
        read_id.to_string(),
        full_read_len,
        anchor.read_start,
        anchor.read_end,
        '+',
        path_dir,
        vec![node_id_num],
        path_offset + r_span,
        path_offset,
        path_offset + r_span,
        match_len,
        match_len.to_string(),
        "60".to_string(),
        comments,
    )
}

fn process_chain_alignment(
    idx: usize,
    read_id_b: &Vec<u8>,
    chain_b: &Chain,
    primary_line: &str,
    reads: &HashMap<String, String>,
    index: &graph_index::PathIndex,
) -> Option<(usize, String, rg_exp::gaf::GAFStruct)> {
    if chain_b.anchors.len() < 2 {
        // Single-anchor chunk: emit minigraph primary alignment as-is in GAF form.
        if chain_b.anchors.len() == 1 {
            let read_id_str = String::from_utf8_lossy(read_id_b).to_string();
            if let Some(gaf) = gaf_from_minigraph_primary_line(primary_line) {
                return Some((idx, read_id_str, gaf));
            }
            eprintln!(
                "Failed to parse single-anchor primary line for read {}: {}",
                String::from_utf8_lossy(read_id_b),
                primary_line
            );
        }
        return None;
    }

    let read_id_str = String::from_utf8_lossy(read_id_b).to_string();
    let read_seq_opt = reads.get(&read_id_str);
    let full_read_len = read_seq_opt.map(|s| s.len()).unwrap_or(0);
    let mut chain_gaf: Option<rg_exp::gaf::GAFStruct> = None;
    // running_path_offset tracks cumulative reference consumed, giving all pieces
    // a consistent global coordinate system so merge_gafs min(path_start) is correct.
    let mut running_path_offset: usize = 0;
    let anchors = &chain_b.anchors;

    for i in 0..anchors.len() {
        let anchor = &anchors[i];

        // Synthetic perfect-match GAF for this anchor (trusted by minigraph).
        let anchor_ref_span = anchor.graph_pos.position.1.saturating_sub(anchor.graph_pos.position.0);
        let anchor_gaf = make_anchor_gaf(anchor, &read_id_str, full_read_len, running_path_offset);
        running_path_offset += anchor_ref_span;
        chain_gaf = Some(match chain_gaf.take() {
            Some(existing) => rg_exp::gaf::GAFStruct::merge_gafs(&existing, &anchor_gaf),
            None => anchor_gaf,
        });

        if i + 1 >= anchors.len() {
            break;
        }
        let from = anchor;
        let to = &anchors[i + 1];

        // align only the gap between the two anchors
        if from.read_end >= to.read_start {
            continue;
        }
        let gap_len = to.read_start - from.read_end;
        if gap_len < 8 {
            continue;
        }

        let mut rc = false;
        let mut read_slice: String;
        if let Some(read_seq) = read_seq_opt {
            read_slice = read_seq[from.read_end..to.read_start].to_string();
            if !from.graph_pos.orientation {
                rc = true;
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
            eprintln!("Read ID {} not found in reads", String::from_utf8_lossy(read_id_b));
            read_slice = String::new();
        }

        let subgraph = chain::extract_subgraph(from, to, index, read_slice.len() + (read_slice.len() / 100));
        let mut gfa_out = String::new();
        gfa::writer::write_gfa(&subgraph, &mut gfa_out);
        gfa_out.insert_str(0, "GRAPH:\n");
        let read_out = format!("READ:\n>{}\n{}", String::from_utf8_lossy(read_id_b), read_slice);
        eprintln!(
            "SUBGRAPH:\n{}\nREAD SLICE ({} to {}):\n>{}\n{}",
            gfa_out, from.read_end, to.read_start,
            String::from_utf8_lossy(read_id_b), read_slice
        );

        let recalign_t0 = std::time::Instant::now();
        let run_recalign = Command::new(RECALIGN_BIN)
            .arg("-efast")
            .arg("-s8")
            .arg("-k2")
            .arg("-r1")
            .arg("-t1")
            .arg("-g")
            .arg("-")
            .arg("-q")
            .arg("-")
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::null())
            .spawn()
            .and_then(|mut child| {
                {
                    let stdin = child.stdin.as_mut().expect("Failed to open stdin");
                    use std::io::Write;
                    stdin.write_all(gfa_out.as_bytes()).expect("Failed to write graph to stdin");
                    stdin.write_all(read_out.as_bytes()).expect("Failed to write read to stdin");
                }
                child.wait_with_output()
            });

        let elapsed = recalign_t0.elapsed().as_secs();
        if elapsed > 10 {
            eprintln!(
                "SLOW recalign: read {} anchors ({}->{}) took {}s",
                String::from_utf8_lossy(read_id_b),
                String::from_utf8_lossy(&from.graph_pos.node_id),
                String::from_utf8_lossy(&to.graph_pos.node_id),
                elapsed
            );
        }

        match run_recalign {
            Ok(recalign_output) => {
                if !recalign_output.status.success() {
                    eprintln!(
                        "Recalign failed for anchor pair ({} -> {}) with error {}, skipping.",
                        String::from_utf8_lossy(&from.graph_pos.node_id),
                        String::from_utf8_lossy(&to.graph_pos.node_id),
                        String::from_utf8_lossy(&recalign_output.stderr)
                    );
                    continue;
                }

                let mut gaf = get_gaf_from_recalign_output(&recalign_output.stdout, rc);
                gaf.query_length = full_read_len;
                gaf.query_start += from.read_end;
                gaf.query_end += from.read_end;

                // Replace recalign's local path coords with global running offset.
                // gap_ref_span is path_end - path_start before any offset replacement,
                // which equals reference_consumed by the CIGAR (consistent after RC too).
                let gap_ref_span = gaf.path_end.saturating_sub(gaf.path_start);
                gaf.path_start = running_path_offset;
                gaf.path_end = running_path_offset + gap_ref_span;
                gaf.path_length = gaf.path_end;
                running_path_offset += gap_ref_span;
                // The RC in get_gaf_from_recalign_output already mapped coordinates and
                // reversed path/CIGAR into the original read's orientation. Force strand='+'
                // so merge_gafs does not apply a second reverse_complement that would
                // re-reverse path and CIGAR and duplicate nodes.
                gaf.strand = '+';

                chain_gaf = Some(match chain_gaf.take() {
                    Some(existing) => rg_exp::gaf::GAFStruct::merge_gafs(&existing, &gaf),
                    None => gaf,
                });
            }
            Err(e) => {
                eprintln!(
                    "Failed to execute recalign at '{}': {}. Build the submodule first.",
                    RECALIGN_BIN, e
                );
            }
        }
    }

    chain_gaf.map(|mut gaf| {
        gaf.reconcile_cigar_with_spans();
        gaf.convert_to_minigraph_comments();
        gaf.strand = '+';
        (idx, read_id_str, gaf)
    })
}

fn gaf_from_minigraph_primary_line(line: &str) -> Option<rg_exp::gaf::GAFStruct> {
    let fields: Vec<&str> = line.trim().split('\t').collect();
    if fields.len() < 12 {
        return None;
    }

    let path_field = fields[5];
    let path_dir = path_field.chars().next().unwrap_or('>');
    let path: Vec<usize> = path_field
        .trim_start_matches(|c| c == '>' || c == '<')
        .split(|c| c == '>' || c == '<')
        .filter(|s| !s.is_empty())
        .filter_map(|id| id.parse::<usize>().ok())
        .collect();

    if path.is_empty() {
        return None;
    }

    let comments = if fields.len() > 12 {
        fields[12..].join("\t")
    } else {
        String::new()
    };

    Some(rg_exp::gaf::GAFStruct::build_gaf_struct(
        fields[0].to_string(),
        fields[1].parse::<usize>().ok()?,
        fields[2].parse::<usize>().ok()?,
        fields[3].parse::<usize>().ok()?,
        fields[4].chars().next().unwrap_or('+'),
        path_dir,
        path,
        fields[6].parse::<usize>().ok()?,
        fields[7].parse::<usize>().ok()?,
        fields[8].parse::<usize>().ok()?,
        fields[9].parse::<usize>().ok()?,
        fields[10].to_string(),
        fields[11].to_string(),
        comments,
    ))
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
