/// GAFStruct represents a gaf alignment, with each field ordered as normal gaf field
#[derive(Debug, Clone)]
pub struct GAFStruct {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub path_dir: char,
    pub path: Vec<usize>,
    pub path_length: usize,
    pub path_start: usize,
    pub path_end: usize,
    pub residue_matches_number: usize,
    pub alignment_block_length: String,
    pub mapping_quality: String,
    pub comments: String,
}

#[derive(Debug, Clone, Copy, Default)]
struct CigarStats {
    query_consumed: usize,
    reference_consumed: usize,
    alignment_block: usize,
    residue_matches: usize,
}

fn parse_cigar_ops(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut n: usize = 0;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            n = n * 10 + c.to_digit(10).unwrap_or(0) as usize;
            continue;
        }
        if n == 0 {
            continue;
        }
        if matches!(c, 'M' | 'I' | 'D' | 'N' | 'S' | 'H' | 'P' | '=' | 'X') {
            ops.push((n, c));
        }
        n = 0;
    }
    ops
}

fn cigar_to_string(ops: &[(usize, char)]) -> String {
    let mut out = String::new();
    for (n, op) in ops {
        out.push_str(&n.to_string());
        out.push(*op);
    }
    out
}

fn merge_adjacent_ops(ops: &[(usize, char)]) -> Vec<(usize, char)> {
    let mut merged: Vec<(usize, char)> = Vec::with_capacity(ops.len());
    for (n, op) in ops {
        if *n == 0 {
            continue;
        }
        if let Some((last_n, last_op)) = merged.last_mut() {
            if *last_op == *op {
                *last_n += *n;
                continue;
            }
        }
        merged.push((*n, *op));
    }
    merged
}

fn trim_cigar_query_prefix(cigar: &str, mut trim_q: usize) -> String {
    if trim_q == 0 || cigar.is_empty() {
        return cigar.to_string();
    }

    let ops = parse_cigar_ops(cigar);
    let mut out: Vec<(usize, char)> = Vec::with_capacity(ops.len());

    for (n, op) in ops {
        if trim_q == 0 {
            out.push((n, op));
            continue;
        }

        let consumes_query = matches!(op, 'M' | '=' | 'X' | 'I');
        if !consumes_query {
            // Deletions before the trim boundary belong to the discarded prefix.
            continue;
        }

        if n <= trim_q {
            trim_q -= n;
        } else {
            out.push((n - trim_q, op));
            trim_q = 0;
        }
    }

    cigar_to_string(&merge_adjacent_ops(&out))
}

fn trim_cigar_suffix_to_spans(cigar: &str, target_q: usize, target_r: usize) -> String {
    let mut ops = parse_cigar_ops(cigar);
    let mut stats = cigar_stats(cigar);

    while (stats.query_consumed > target_q || stats.reference_consumed > target_r) && !ops.is_empty() {
        let (n, op) = *ops.last().unwrap();
        let q_cons = matches!(op, 'M' | '=' | 'X' | 'I');
        let r_cons = matches!(op, 'M' | '=' | 'X' | 'D' | 'N');

        let dq = stats.query_consumed.saturating_sub(target_q);
        let dr = stats.reference_consumed.saturating_sub(target_r);

        let removable = match (q_cons, r_cons) {
            (true, true) if dq > 0 && dr > 0 => n.min(dq).min(dr),
            (true, false) if dq > 0 => n.min(dq),
            (false, true) if dr > 0 => n.min(dr),
            (false, false) => n,
            _ => 0,
        };

        if removable == 0 {
            break;
        }

        let mut popped = ops.pop().unwrap();
        popped.0 -= removable;
        if popped.0 > 0 {
            ops.push(popped);
        }

        if q_cons {
            stats.query_consumed = stats.query_consumed.saturating_sub(removable);
        }
        if r_cons {
            stats.reference_consumed = stats.reference_consumed.saturating_sub(removable);
        }
        if matches!(op, 'M' | '=' | 'X' | 'I' | 'D' | 'N') {
            stats.alignment_block = stats.alignment_block.saturating_sub(removable);
        }
        if matches!(op, 'M' | '=') {
            stats.residue_matches = stats.residue_matches.saturating_sub(removable);
        }
    }

    cigar_to_string(&merge_adjacent_ops(&ops))
}

fn adjust_cigar_to_spans(cigar: &str, target_q: usize, target_r: usize) -> String {
    let mut adjusted = trim_cigar_suffix_to_spans(cigar, target_q, target_r);
    let s = cigar_stats(&adjusted);

    let mut tail: Vec<(usize, char)> = Vec::new();
    let dq = target_q.saturating_sub(s.query_consumed);
    let dr = target_r.saturating_sub(s.reference_consumed);

    let m = dq.min(dr);
    if m > 0 {
        tail.push((m, 'M'));
    }
    if dq > dr {
        tail.push((dq - dr, 'I'));
    } else if dr > dq {
        tail.push((dr - dq, 'D'));
    }

    if !tail.is_empty() {
        let mut ops = parse_cigar_ops(&adjusted);
        ops.extend(tail);
        adjusted = cigar_to_string(&merge_adjacent_ops(&ops));
    }

    let mut ops = parse_cigar_ops(&adjusted);
    let mut s2 = cigar_stats(&adjusted);
    if s2.reference_consumed > target_r {
        let mut trim_r = s2.reference_consumed - target_r;
        let mut q_removed = 0usize;

        let mut i = ops.len();
        while trim_r > 0 && i > 0 {
            i -= 1;
            let (n, op) = ops[i];
            let r_cons = matches!(op, 'M' | '=' | 'X' | 'D' | 'N');
            let q_cons = matches!(op, 'M' | '=' | 'X' | 'I');
            if !r_cons || n == 0 {
                continue;
            }

            let take = n.min(trim_r);
            ops[i].0 -= take;
            trim_r -= take;
            if q_cons {
                q_removed += take;
            }
        }

        ops.retain(|(n, _)| *n > 0);
        if q_removed > 0 {
            ops.push((q_removed, 'I'));
        }
        adjusted = cigar_to_string(&merge_adjacent_ops(&ops));
        s2 = cigar_stats(&adjusted);
    }

    if s2.query_consumed < target_q {
        let mut ops2 = parse_cigar_ops(&adjusted);
        ops2.push((target_q - s2.query_consumed, 'I'));
        adjusted = cigar_to_string(&merge_adjacent_ops(&ops2));
    }

    adjusted
}

fn reverse_cigar(cigar: &str) -> String {
    let mut ops = parse_cigar_ops(cigar);
    ops.reverse();
    cigar_to_string(&ops)
}

fn cigar_stats(cigar: &str) -> CigarStats {
    let mut s = CigarStats::default();
    for (n, op) in parse_cigar_ops(cigar) {
        match op {
            'M' => {
                s.query_consumed += n;
                s.reference_consumed += n;
                s.alignment_block += n;
                s.residue_matches += n;
            }
            '=' => {
                s.query_consumed += n;
                s.reference_consumed += n;
                s.alignment_block += n;
                s.residue_matches += n;
            }
            'X' => {
                s.query_consumed += n;
                s.reference_consumed += n;
                s.alignment_block += n;
            }
            'I' => {
                s.query_consumed += n;
                s.alignment_block += n;
            }
            'D' | 'N' => {
                s.reference_consumed += n;
                s.alignment_block += n;
            }
            _ => {}
        }
    }
    s
}

fn nm_from_cigar(cigar: &str) -> usize {
    parse_cigar_ops(cigar)
        .iter()
        .filter_map(|(n, op)| if matches!(op, 'X' | 'I' | 'D') { Some(*n) } else { None })
        .sum()
}

fn extract_cigar_from_comments(comments: &str) -> Option<String> {
    if let Some(tag) = comments.split('\t').find(|t| t.starts_with("cg:Z:")) {
        return Some(tag[5..].to_string());
    }
    if comments.starts_with("@CO,") {
        let parts: Vec<&str> = comments.split(',').collect();
        if parts.len() >= 3 {
            return Some(parts[2].to_string());
        }
    }
    None
}

fn replace_cigar_in_comments(comments: &str, new_cigar: &str) -> String {
    if comments.starts_with("@CO,") {
        let mut parts: Vec<String> = comments.split(',').map(|s| s.to_string()).collect();
        if parts.len() >= 3 {
            parts[2] = new_cigar.to_string();
            return parts.join(",");
        }
    }

    let mut tags: Vec<String> = if comments.is_empty() {
        Vec::new()
    } else {
        comments.split('\t').map(|s| s.to_string()).collect()
    };

    if let Some(i) = tags.iter().position(|t| t.starts_with("cg:Z:")) {
        tags[i] = format!("cg:Z:{}", new_cigar);
    } else {
        tags.push(format!("cg:Z:{}", new_cigar));
    }
    tags.join("\t")
}

impl Default for GAFStruct {
    fn default() -> Self {
        Self::new()
    }
}

impl GAFStruct {
    pub fn new() -> GAFStruct {
        GAFStruct {
            query_name: String::from(""),
            query_length: 0,
            query_start: 0,
            query_end: 0,
            strand: ' ',
            path_dir: '>',
            path: vec![0usize],
            path_length: 0,
            path_start: 0,
            path_end: 0,
            residue_matches_number: 0,
            alignment_block_length: String::from(""),
            mapping_quality: String::from(""),
            comments: String::from(""),
        }
    }

    pub fn from_string(gaf_line: &str) -> GAFStruct {
        let fields: Vec<&str> = gaf_line.trim().split('\t').collect();
        let path_dir = fields[5].chars().next().unwrap_or('>');
        let path_ids: Vec<usize> = fields[5]
            .trim_start_matches(|c| c == '>' || c == '<')
            .split(|c| c == '>' || c == '<')
            .filter(|s| !s.is_empty())
            .map(|id| id.parse::<usize>().unwrap())
            .collect();
        GAFStruct {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse::<usize>().unwrap(),
            query_start: fields[2].parse::<usize>().unwrap(),
            query_end: fields[3].parse::<usize>().unwrap(),
            strand: fields[4].chars().next().unwrap(),
            path_dir,
            path: path_ids,
            path_length: fields[6].parse::<usize>().unwrap(),
            path_start: fields[7].parse::<usize>().unwrap(),
            path_end: fields[8].parse::<usize>().unwrap(),
            residue_matches_number: fields[9].parse::<usize>().unwrap(),
            alignment_block_length: fields[10].to_string(),
            mapping_quality: fields[11].to_string(),
            comments: fields[12].to_string(),
        }
    }

    pub fn reverse_complement(&self) -> GAFStruct {
        let new_strand = if self.strand == '+' { '-' } else { '+' };
        let new_path_dir = if self.path_dir == '>' { '<' } else { '>' };
        let new_path: Vec<usize> = self.path.iter().rev().cloned().collect();

        let new_comments = if let Some(cigar) = extract_cigar_from_comments(&self.comments) {
            replace_cigar_in_comments(&self.comments, &reverse_cigar(&cigar))
        } else {
            self.comments.clone()
        };

        GAFStruct {
            query_name: self.query_name.clone(),
            query_length: self.query_length,
            query_start: self.query_length - self.query_end,
            query_end: self.query_length - self.query_start,
            strand: new_strand,
            path_dir: new_path_dir,
            path: new_path,
            path_length: self.path_length,
            path_start: self.path_length - self.path_end,
            path_end: self.path_length - self.path_start,
            residue_matches_number: self.residue_matches_number,
            alignment_block_length: self.alignment_block_length.clone(),
            mapping_quality: self.mapping_quality.clone(),
            comments: new_comments,
        }
    }
    pub fn merge_gafs(gaf1: &GAFStruct, gaf2: &GAFStruct) -> GAFStruct {
        let gaf2_adj = if gaf1.strand != gaf2.strand {
            gaf2.reverse_complement()
        } else {
            gaf2.clone()
        };

        let mut new_path: Vec<usize> = Vec::with_capacity(gaf1.path.len() + gaf2_adj.path.len());
        new_path.extend_from_slice(&gaf1.path);
        if !gaf2_adj.path.is_empty() {
            if gaf1.path.last() == gaf2_adj.path.first() {
                // contiguous on nodes: skip first node of gaf2_adj
                new_path.extend_from_slice(&gaf2_adj.path[1..]);
            } else {
                new_path.extend_from_slice(&gaf2_adj.path);
            }
        }

        let new_query_start = std::cmp::min(gaf1.query_start, gaf2_adj.query_start);
        let new_query_end = std::cmp::max(gaf1.query_end, gaf2_adj.query_end);

        let c1 = extract_cigar_from_comments(&gaf1.comments).unwrap_or_default();
        let c2 = extract_cigar_from_comments(&gaf2_adj.comments).unwrap_or_default();

        let overlap_q = gaf1.query_end.saturating_sub(gaf2_adj.query_start);
        let c2_trimmed = trim_cigar_query_prefix(&c2, overlap_q);

        let mut merged_ops = parse_cigar_ops(&c1);
        merged_ops.extend(parse_cigar_ops(&c2_trimmed));
        let merged_cigar = cigar_to_string(&merge_adjacent_ops(&merged_ops));
        let merged_comments = if merged_cigar.is_empty() {
            if !gaf2_adj.comments.is_empty() {
                gaf2_adj.comments.clone()
            } else {
                gaf1.comments.clone()
            }
        } else {
            replace_cigar_in_comments(&gaf1.comments, &merged_cigar)
        };

        let stats = cigar_stats(&merged_cigar);

        let new_path_start = std::cmp::min(gaf1.path_start, gaf2_adj.path_start);
        let new_path_end = new_path_start + stats.reference_consumed;
        let new_path_length = std::cmp::max(
            std::cmp::max(gaf1.path_length, gaf2_adj.path_length),
            new_path_end,
        );

        let new_residue_matches = if stats.residue_matches > 0 {
            stats.residue_matches
        } else {
            gaf1.residue_matches_number + gaf2_adj.residue_matches_number
        };
        let new_alignment_block_length = if stats.alignment_block > 0 {
            stats.alignment_block.to_string()
        } else {
            (gaf1.alignment_block_length.parse::<usize>().unwrap_or(0)
                + gaf2_adj.alignment_block_length.parse::<usize>().unwrap_or(0))
                .to_string()
        };

        let new_mapping_quality = match (gaf1.mapping_quality.parse::<usize>(), gaf2_adj.mapping_quality.parse::<usize>()) {
            (Ok(a), Ok(b)) => std::cmp::max(a, b).to_string(),
            _ => {
                if !gaf1.mapping_quality.is_empty() {
                    gaf1.mapping_quality.clone()
                } else {
                    gaf2_adj.mapping_quality.clone()
                }
            }
        };

        GAFStruct::build_gaf_struct(
            gaf1.query_name.clone(),
            gaf1.query_length, // query_length is the read length, invariant across merges
            new_query_start,
            new_query_end,
            gaf1.strand,
            gaf1.path_dir,
            new_path,
            new_path_length,
            new_path_start,
            new_path_end,
            new_residue_matches,
            new_alignment_block_length,
            new_mapping_quality,
            merged_comments,
        )
    }

    pub fn reconcile_cigar_with_spans(&mut self) {
        let cigar = match extract_cigar_from_comments(&self.comments) {
            Some(c) => c,
            None => return,
        };

        let target_q = self.query_end.saturating_sub(self.query_start);
        let target_r = self.path_end.saturating_sub(self.path_start);
        let adjusted = adjust_cigar_to_spans(&cigar, target_q, target_r);
        self.comments = replace_cigar_in_comments(&self.comments, &adjusted);

        let stats = cigar_stats(&adjusted);
        self.residue_matches_number = stats.residue_matches;
        self.alignment_block_length = stats.alignment_block.to_string();
    }

    pub fn build_gaf_struct(
        query_name: String,
        query_length: usize,
        query_start: usize,
        query_end: usize,
        strand: char,
        path_dir: char,
        path: Vec<usize>,
        path_length: usize,
        path_start: usize,
        path_end: usize,
        residue_matches_number: usize,
        alignment_block_length: String,
        mapping_quality: String,
        comments: String,
    ) -> GAFStruct {
        GAFStruct {
            query_name,
            query_length,
            query_start,
            query_end,
            strand,
            path_dir,
            path,
            path_length,
            path_start,
            path_end,
            residue_matches_number,
            alignment_block_length,
            mapping_quality,
            comments,
        }
    }
    // Converts recalign's @CO,... comment format to standard minigraph-style SAM tags.
    // No-ops if comments are already in tag format (e.g. from a minigraph primary line).
    pub fn convert_to_minigraph_comments(&mut self) {
        if !self.comments.starts_with("@CO,") {
            return;
        }
        let raw_cigar = extract_cigar_from_comments(&self.comments).unwrap_or_default();
        let cigar = raw_cigar.replace('M', "="); // recalign M == exact match
        let nm = nm_from_cigar(&cigar);
        let block = self.alignment_block_length.parse::<usize>().unwrap_or(0);
        let dv = if block > 0 { nm as f64 / block as f64 } else { 0.0 };
        let mut tags = format!("tp:A:P\tNM:i:{}\tdv:f:{:.4}", nm, dv);
        if !cigar.is_empty() {
            tags.push_str(&format!("\tcg:Z:{}", cigar));
        }
        self.comments = tags;
    }

    pub fn to_string(self) -> String {
        let dir = if self.path_dir == '>' { ">" } else { "<" };
        let path_matching: String = self
            .path
            .iter()
            .map(|id| id.to_string())
            .collect::<Vec<String>>()
            .join(dir);
        let gaf_struct_to_string = format!(
            "{}\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            dir,
            path_matching,
            self.path_length,
            self.path_start,
            self.path_end,
            self.residue_matches_number,
            self.alignment_block_length,
            self.mapping_quality,
            self.comments
        );
        gaf_struct_to_string
    }
}