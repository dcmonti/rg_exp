use gfa::gfa::GFA;
use std::{collections::HashMap, str};

pub fn index_paths(gfa_file: &str) -> PathIndex {
    let parser = gfa::parser::GFAParser::new();
    let gfa: GFA<Vec<u8>, ()> = if gfa_file.ends_with(".gz") {
        use flate2::read::MultiGzDecoder;
        use std::{fs::File, io::Read};
        let f = File::open(gfa_file).expect("Cannot open GFA file");
        let mut buf = Vec::new();
        MultiGzDecoder::new(f)
            .read_to_end(&mut buf)
            .expect("Failed to decompress GFA");
        parser.parse_lines(buf.split(|&b| b == b'\n')).unwrap()
    } else {
        parser.parse_file(gfa_file).unwrap()
    };
    let mut index = PathIndex::new();
    for segment in &gfa.segments {
        index.add_node(segment.name.clone(), segment.sequence.clone());
    }

    for p in &gfa.paths {
        let new_path_id = index.paths_to_ids.len();
        index.paths_to_ids.insert(new_path_id, p.path_name.clone());
        let mut nodes_in_path = Vec::new();
        p.iter()
            .enumerate()
            .for_each(|(position, (node_id, orientation))| {
                index.update_path_occurrence(node_id, new_path_id, position);
                let dir = if orientation == gfa::gfa::Orientation::Forward {
                    true
                } else {
                    false
                };
                nodes_in_path.push((node_id.to_vec(), dir));
            });
        index.add_path(new_path_id, p.path_name.clone(), nodes_in_path);
    }
    eprintln!("Indexed {} paths.", index.paths.len());
    index
}

pub struct NodeInfo {
    pub segment: Vec<u8>,
    pub path_positions: HashMap<usize, (usize, usize)>, // path_id -> (first_occ, last_occ)
}

impl NodeInfo {
    fn new(segment: Vec<u8>) -> Self {
        NodeInfo {
            segment,
            path_positions: HashMap::new(),
        }
    }

    fn add_path_occurrence(&mut self, path_id: usize, position: usize) {
        self.path_positions
            .entry(path_id)
            .and_modify(|pos| pos.1 = position)
            .or_insert((position, position));
    }

    pub fn get_path_fl(&self, path_id: &usize) -> Option<&(usize, usize)> {
        self.path_positions.get(path_id)
    }

    pub fn get_paths(&self) -> Vec<usize> {
        self.path_positions.keys().cloned().collect()
    }
}
pub struct PathIndex {
    pub nodes: HashMap<Vec<u8>, NodeInfo>,
    pub paths: HashMap<usize, PathInfo>,
    pub paths_to_ids: HashMap<usize, Vec<u8>>,
}

pub struct PathInfo {
    pub path_name: Vec<u8>,
    pub nodes: Vec<(Vec<u8>, bool)>,
}
impl PathInfo {
    pub fn new(path_name: Vec<u8>, nodes: Vec<(Vec<u8>, bool)>) -> Self {
        PathInfo { path_name, nodes }
    }

    pub fn iter(&self) -> impl Iterator<Item = &(Vec<u8>, bool)> {
        self.nodes.iter()
    }
}

pub enum PathRelation {
    First,
    Second,
    Both,
}

impl PathIndex {
    fn new() -> Self {
        PathIndex {
            nodes: HashMap::new(),
            paths: HashMap::new(),
            paths_to_ids: HashMap::new(),
        }
    }

    fn add_node(&mut self, node_id: Vec<u8>, segment: Vec<u8>) {
        self.nodes.insert(node_id, NodeInfo::new(segment));
    }

    fn update_path_occurrence(&mut self, node_id: &[u8], path_id: usize, position: usize) {
        if let Some(node_info) = self.nodes.get_mut(node_id) {
            node_info.add_path_occurrence(path_id, position);
        }
    }
    fn add_path(&mut self, path_id: usize, path_name: Vec<u8>, nodes: Vec<(Vec<u8>, bool)>) {
        self.paths.insert(path_id, PathInfo { path_name, nodes });
    }

    fn get_node_path_fl(&self, node_id: &[u8], path_id: &usize) -> Option<&(usize, usize)> {
        self.nodes
            .get(node_id)
            .and_then(|node_info| node_info.get_path_fl(path_id))
    }

    pub fn get(&self, node_id: &[u8]) -> Option<&NodeInfo> {
        self.nodes.get(node_id)
    }

    pub fn common_paths(&self, node_a: &[u8], node_b: &[u8]) -> Vec<(usize, PathRelation)> {
        let paths_a = match self.get(node_a) {
            Some(info) => info.get_paths(),
            None => return Vec::new(),
        };
        let paths_b = match self.get(node_b) {
            Some(info) => info.get_paths(),
            None => return Vec::new(),
        };
        paths_a
            .into_iter()
            .filter(|path_id| paths_b.contains(path_id))
            .map(|p| (p, PathRelation::Both))
            .collect()
        
    }

    pub fn both_paths(&self, node_a: &[u8], node_b: &[u8]) -> Vec<(usize, PathRelation)> {
        let paths_a = match self.get(node_a) {
            Some(info) => info.get_paths(),
            None => return Vec::new(),
        };
        let paths_b = match self.get(node_b) {
            Some(info) => info.get_paths(),
            None => return Vec::new(),
        };
        let mut both_paths = paths_a.clone();
        both_paths.extend(paths_b.clone());
        both_paths.sort_unstable();
        both_paths.dedup();
        // set if paths are in both nodes or only in one
        let annotated_paths = both_paths
            .into_iter()
            .map(|p| {if paths_a.contains(&p) && paths_b.contains(&p) {
                (p, PathRelation::Both)
            } else if paths_a.contains(&p) {
                (p, PathRelation::First)
            } else {
                (p, PathRelation::Second)
            }})
            .collect::<Vec<(usize, PathRelation)>>();
        annotated_paths
    }
    fn serialize_index(&self) -> Vec<u8> {
        let mut serialized = Vec::new();
        for (node_id, node_info) in &self.nodes {
            serialized.extend_from_slice(node_id);
            for (path_id, (first_occ, last_occ)) in &node_info.path_positions {
                serialized.extend_from_slice(&path_id.to_le_bytes());
                serialized.extend_from_slice(&first_occ.to_le_bytes());
                serialized.extend_from_slice(&last_occ.to_le_bytes());
            }
        }
        serialized
    }

    pub fn dump_index(&self, path: &str) {
        let serialized = self.serialize_index();
        std::fs::write(path, serialized).expect("Unable to write index to file");
    }

    fn deserialize_index(data: &[u8]) -> Self {
        let mut index = PathIndex::new();
        let mut i = 0;
        while i < data.len() {
            let node_id = data[i..i + 16].to_vec(); // assuming fixed length for simplicity
            i += 16;
            let mut node_info = NodeInfo::new(Vec::new());
            while i + 24 <= data.len() {
                let path_id = usize::from_le_bytes(data[i..i + 8].try_into().unwrap());
                let first_occ = usize::from_le_bytes(data[i + 8..i + 16].try_into().unwrap());
                let last_occ = usize::from_le_bytes(data[i + 16..i + 24].try_into().unwrap());
                node_info
                    .path_positions
                    .insert(path_id, (first_occ, last_occ));
                i += 24;
            }
            index.nodes.insert(node_id, node_info);
        }
        index
    }

    pub fn load_index(path: &str) -> Self {
        let data = std::fs::read(path).expect("Unable to read index from file");
        Self::deserialize_index(&data)
    }
}
