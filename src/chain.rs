use std::fmt::Display;

use gfa::gfa::GFA;
use crate::graph_index;

pub fn build_chain(alignment: &str) -> Chain {
    // Chain from minigraph alignment output parsing logic goes here
    // Expected input:
    // * sName sLen nMz div sStart sEnd qStart qEnd\n
    let mut chain = Chain::new(vec![]);
    alignment.lines().for_each(|line| {
        let fields = line.split('\t').collect::<Vec<&str>>();
        if fields.len() == 9 {
            let s_name = &fields[1][1..];
            let s_start: usize = fields[5].parse().unwrap();
            let s_end: usize = fields[6].parse().unwrap();
            let q_start: usize = fields[7].parse().unwrap();
            let q_end: usize = fields[8].parse().unwrap();
            let orientation = fields[1].chars().next().unwrap() == '>';
            let graph_pos = Pos::new(s_name.as_bytes().to_vec(), (s_start, s_end), orientation);
            let anchor = Anchor::new(graph_pos, q_start, q_end);
            chain.add_anchor(anchor);
        }
    });
    chain
}

pub fn extract_subgraph(
    from_anchor: &Anchor,
    to_anchor: &Anchor,
    path_index: &graph_index::PathIndex,
    substr_len: usize,
) -> GFA<Vec<u8>, ()> {
    let from_node = &from_anchor.graph_pos.node_id;
    let to_node = &to_anchor.graph_pos.node_id;
    // paths to keep are those either in from or to node
    let paths_to_keep = path_index.both_paths(from_node, to_node);
    // keep only segments and links that are part of the paths in path_names_to_keep between from_node and to_node
    let mut subgraph: GFA<Vec<u8>, ()> = GFA::new();
    let mut segments_set = std::collections::HashSet::new();
    // Further implementation needed to build the subgraph
    for (path_id, relation) in paths_to_keep {
        let path_name = &path_index.paths_to_ids[&path_id];
        
        // compute start and end indices in the path, 
        // if in Both relation, use from_node and to_node positions
        // if in First relation, use from_node position to first node at length substr_len
        // if in Second relation, use to_node position back to first node at length substr_len
        let (st_idx, en_idx) = match relation {
            graph_index::PathRelation::Both => {
                let st = path_index
                    .get(&from_node)
                    .unwrap()
                    .get_path_fl(&path_id)
                    .unwrap()
                    .0;
                let en = path_index
                    .get(&to_node)
                    .unwrap()
                    .get_path_fl(&path_id)
                    .unwrap()
                    .1;
                (st, en)
            },
            graph_index::PathRelation::First => {
                let st = path_index
                    .get(&from_node)
                    .unwrap()
                    .get_path_fl(&path_id)
                    .unwrap()
                    .0;
                let mut en = st;
                let mut curr_len = 0;
                while curr_len < substr_len && en < path_index.paths[&path_id].nodes.len() {
                    let next_node = &path_index.paths[&path_id].nodes[en].0;
                    let next_node_len = path_index.get(&next_node).unwrap().segment.len();
                    curr_len += next_node_len;
                    en += 1;
                }
                (st, en)
            },
            graph_index::PathRelation::Second => {
                let en = path_index
                    .get(&to_node)
                    .unwrap()
                    .get_path_fl(&path_id)
                    .unwrap()
                    .1;
                let mut st = en;
                let mut curr_len = 0;
                while curr_len < substr_len && st > 0 {
                    let prev_node = &path_index.paths[&path_id].nodes[st].0;
                    let prev_node_len = path_index.get(&prev_node).unwrap().segment.len();
                    curr_len += prev_node_len;
                    if st == 0 {
                        break;
                    }
                    st -= 1;
                }
                (st, en)
            },
        };

        // Add segments between from_node and to_node in this path
        let mut path_nodes = Vec::new();
        let nodes = &path_index.paths[&path_id].nodes;
        // ensure indices are within bounds
        let mut st_idx = std::cmp::min(st_idx, nodes.len().saturating_sub(1));
        let mut en_idx = std::cmp::min(en_idx, nodes.len().saturating_sub(1));
        if st_idx > en_idx {
           let tmp = st_idx;
           st_idx = en_idx;
           en_idx = tmp;
        }
        nodes[st_idx..=en_idx].iter().for_each(|(node_id, orientation)| {
            let mut seq_slice: &[u8] = path_index.get(&node_id).unwrap().segment.as_slice();
            // if first node, trim to from position start
            if node_id == from_node {
                let from_pos = from_anchor.graph_pos.position.0;
                seq_slice = &seq_slice[from_pos..];
            }
            // if last node, trim to to position end
            if node_id == to_node {
                let to_pos = to_anchor.graph_pos.position.1;
                seq_slice = &seq_slice[..to_pos];
            }
            let id_slice: &[u8] = node_id.as_slice();
            let seg = gfa::gfa::Segment {
                name: id_slice.to_vec(),
                sequence: seq_slice.to_vec(),
                ..Default::default()
            };
            if !segments_set.contains(id_slice) {
                segments_set.insert(id_slice.to_vec());
                subgraph.segments.push(seg);
            }
            let direction = if *orientation {
                [b'+', b',']
            } else {
                [b'-', b',']
            };
            let mut new_node = node_id.to_vec();
            new_node.extend_from_slice(&direction);
            path_nodes.extend_from_slice(&new_node);

        }); 
        path_nodes.pop(); // remove trailing comma
        //add \t* at the end
        path_nodes.extend_from_slice(b"\t*");
        let new_path: gfa::gfa::Path<Vec<u8>, _> = gfa::gfa::Path::new(
            path_name.clone(),
            path_nodes,
            Vec::new(),
            Default::default(),
        );
        subgraph.paths.push(new_path);
    }
    subgraph
}

#[derive(Debug)]
pub struct Pos {
    pub node_id: Vec<u8>,
    pub position: (usize, usize),
    pub orientation: bool,
}

impl Pos {
    pub fn new(node_id: Vec<u8>, position: (usize, usize), orientation: bool) -> Self {
        Pos {
            node_id,
            position,
            orientation,
        }
    }
}

impl Display for Pos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Pos(node_id: {}, position: ({}, {}), orientation: {})",
            String::from_utf8_lossy(&self.node_id), self.position.0, self.position.1, self.orientation
        )
    }
}

#[derive(Debug)]
pub struct Anchor {
    pub graph_pos: Pos,
    pub read_start: usize,
    pub read_end: usize,
}

impl Anchor {
    pub fn new(graph_pos: Pos, read_start: usize, read_end: usize) -> Self {
        Anchor {
            graph_pos,
            read_start,
            read_end,
        }
    }
}


impl Display for Anchor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Anchor(graph_pos: {}, read_start: {}, read_end: {})",
            self.graph_pos, self.read_start, self.read_end
        )
    }
}
#[derive(Debug)]
pub struct Chain {
    pub anchors: Vec<Anchor>,
}

impl Chain {
    pub fn new(anchors: Vec<Anchor>) -> Self {
        Chain { anchors }
    }

    pub fn add_anchor(&mut self, anchor: Anchor) {
        self.anchors.push(anchor);
    }
}
