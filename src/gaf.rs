/// GAFStruct represents a gaf alignment, with each field ordered as normal gaf field
#[derive(Debug, Clone)]
pub struct GAFStruct {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub path: Vec<usize>,
    pub path_length: usize,
    pub path_start: usize,
    pub path_end: usize,
    pub residue_matches_number: usize,
    pub alignment_block_length: String,
    pub mapping_quality: String,
    pub comments: String,
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
        let path_ids: Vec<usize> = fields[5]
            .trim_start_matches('>')
            .split('>')
            .map(|id| id.parse::<usize>().unwrap())
            .collect();
        GAFStruct {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse::<usize>().unwrap(),
            query_start: fields[2].parse::<usize>().unwrap(),
            query_end: fields[3].parse::<usize>().unwrap(),
            strand: fields[4].chars().next().unwrap(),
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
        // flip path 
        let new_path: Vec<usize> = self.path.iter().rev().cloned().collect();
        GAFStruct {
            query_name: self.query_name.clone(),
            query_length: self.query_length,
            query_start: self.query_start,
            query_end: self.query_end,
            strand: new_strand,
            path: new_path,
            path_length: self.path_length,
            path_start: self.path_start,
            path_end: self.path_end,
            residue_matches_number: self.residue_matches_number,
            alignment_block_length: self.alignment_block_length.clone(),
            mapping_quality: self.mapping_quality.clone(),
            comments: self.comments.clone(),
        }
    }
    pub fn merge_gafs(gaf1: &GAFStruct, gaf2: &GAFStruct) -> GAFStruct {
        let gaf2_adj = if gaf1.strand != gaf2.strand {
            gaf2.reverse_complement()
        } else {
            gaf2.clone()
        };

        let mut meta_comments = String::new();
        if gaf1.query_name != gaf2_adj.query_name {
            meta_comments.push_str(&format!("MERGE_WARNING: different query_name ({} vs {}); ", gaf1.query_name, gaf2_adj.query_name));
        }
        let contiguity = if gaf1.query_end == gaf2_adj.query_start {
            "contiguous"
        } else if gaf1.query_end > gaf2_adj.query_start {
            meta_comments.push_str(&format!("MERGE_INFO: overlapping by {} bases; ", gaf1.query_end - gaf2_adj.query_start));
            "overlap"
        } else {
            meta_comments.push_str(&format!("MERGE_INFO: gap of {} bases; ", gaf2_adj.query_start - gaf1.query_end));
            "gap"
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
        let new_path_start = std::cmp::min(gaf1.path_start, gaf2_adj.path_start);
        let new_path_end = std::cmp::max(gaf1.path_end, gaf2_adj.path_end);

        let new_residue_matches = gaf1.residue_matches_number + gaf2_adj.residue_matches_number;
        if contiguity != "contiguous" {
            meta_comments.push_str("MERGE_NOTE: residue_matches were summed without overlap-disambiguation; ");
        }
        let new_block_length = gaf1.alignment_block_length.parse::<usize>().unwrap_or(0)
            + gaf2_adj.alignment_block_length.parse::<usize>().unwrap_or(0);
        let new_alignment_block_length = format!("{}", new_block_length);

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

        let new_comments = format!(
            "MERGE_FRAMEWORK; contiguity={}; gaf1_comments={:?}; gaf2_comments={:?}; {}",
            contiguity, gaf1.comments, gaf2_adj.comments, meta_comments
        );

        GAFStruct::build_gaf_struct(
            gaf1.query_name.clone(),
            gaf1.query_length + gaf2_adj.query_length,
            new_query_start,
            new_query_end,
            gaf1.strand,
            new_path,
            gaf1.path_length + gaf2_adj.path_length,
            new_path_start,
            new_path_end,
            new_residue_matches,
            new_alignment_block_length,
            new_mapping_quality,
            //new_comments,
            gaf2_adj.comments.clone()
        )
    }

    pub fn build_gaf_struct(
        query_name: String,
        query_length: usize,
        query_start: usize,
        query_end: usize,
        strand: char,
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
    pub fn to_string(self) -> String {
        let dir = if self.strand == '+' { ">" } else { "<" };
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
