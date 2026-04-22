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
    // reversed = anchors map to the reverse strand; from_node comes AFTER to_node in path order
    let reversed = !from_anchor.graph_pos.orientation;
    let paths_to_keep = path_index.both_paths(from_node, to_node);
    let mut subgraph: GFA<Vec<u8>, ()> = GFA::new();
    let mut segments_set = std::collections::HashSet::new();

    for (path_id, relation) in paths_to_keep {
        let path_name = &path_index.paths_to_ids[&path_id];
        let path_nodes_list = &path_index.paths[&path_id].nodes;

        // Compute (st_idx, en_idx) as a sorted [lo, hi] range to slice path_nodes_list.
        // `reversed` controls iteration order and orientation flipping later.
        let (st_idx, en_idx) = match relation {
            graph_index::PathRelation::Both => {
                let from_idx = path_index.get(from_node).unwrap().get_path_fl(&path_id).unwrap().0;
                let to_idx   = path_index.get(to_node).unwrap().get_path_fl(&path_id).unwrap().1;
                // For reversed reads from_node is after to_node in path order, so from_idx > to_idx.
                // Normalise to lo..=hi; iteration direction is decided below.
                if from_idx <= to_idx { (from_idx, to_idx) } else { (to_idx, from_idx) }
            },
            graph_index::PathRelation::First => {
                // Only from_node is in this path.
                let anchor_idx = path_index.get(from_node).unwrap().get_path_fl(&path_id).unwrap().0;
                if reversed {
                    // Reversed: walk backward (toward lower indices) from from_node.
                    let mut st = anchor_idx;
                    let mut curr_len = 0;
                    while curr_len < substr_len && st > 0 {
                        curr_len += path_index.get(&path_nodes_list[st].0).unwrap().segment.len();
                        st -= 1;
                    }
                    (st, anchor_idx)
                } else {
                    // Forward: walk forward (toward higher indices) from from_node.
                    let mut en = anchor_idx;
                    let mut curr_len = 0;
                    while curr_len < substr_len && en < path_nodes_list.len() {
                        curr_len += path_index.get(&path_nodes_list[en].0).unwrap().segment.len();
                        en += 1;
                    }
                    (anchor_idx, en)
                }
            },
            graph_index::PathRelation::Second => {
                // Only to_node is in this path.
                let anchor_idx = path_index.get(to_node).unwrap().get_path_fl(&path_id).unwrap().1;
                if reversed {
                    // Reversed: walk forward (toward higher indices) from to_node.
                    let mut en = anchor_idx;
                    let mut curr_len = 0;
                    while curr_len < substr_len && en < path_nodes_list.len() {
                        curr_len += path_index.get(&path_nodes_list[en].0).unwrap().segment.len();
                        en += 1;
                    }
                    (anchor_idx, en)
                } else {
                    // Forward: walk backward (toward lower indices) to to_node.
                    let mut st = anchor_idx;
                    let mut curr_len = 0;
                    while curr_len < substr_len && st > 0 {
                        curr_len += path_index.get(&path_nodes_list[st].0).unwrap().segment.len();
                        st -= 1;
                    }
                    (st, anchor_idx)
                }
            },
        };

        let st_idx = st_idx.min(path_nodes_list.len().saturating_sub(1));
        let en_idx = en_idx.min(path_nodes_list.len().saturating_sub(1));
        let (st_idx, en_idx) = if st_idx <= en_idx { (st_idx, en_idx) } else { (en_idx, st_idx) };

        let mut path_nodes = Vec::new();

        // Always iterate in path order. For reversed reads, from_node sits at the high
        // end and to_node at the low end of the range; the RC'd read (sent to recalign)
        // aligns left-to-right against the forward-orientation path segments.
        for (node_id, path_orientation) in &path_nodes_list[st_idx..=en_idx] {
            let mut seq_slice: &[u8] = path_index.get(node_id).unwrap().segment.as_slice();

            if node_id == from_node {
                if reversed {
                    // from_node is the LAST segment for the RC'd read; trim like a to_node
                    // in the forward case: keep [..s_end] so the gap toward the next anchor
                    // is included.
                    seq_slice = &seq_slice[..from_anchor.graph_pos.position.1];
                } else {
                    seq_slice = &seq_slice[from_anchor.graph_pos.position.0..];
                }
            }
            if node_id == to_node {
                if reversed {
                    // to_node is the FIRST segment for the RC'd read; trim like a from_node
                    // in the forward case: keep [s_start..] so the gap is included.
                    seq_slice = &seq_slice[to_anchor.graph_pos.position.0..];
                } else {
                    seq_slice = &seq_slice[..to_anchor.graph_pos.position.1];
                }
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

            let direction = if *path_orientation { [b'+', b','] } else { [b'-', b','] };
            let mut new_node = node_id.to_vec();
            new_node.extend_from_slice(&direction);
            path_nodes.extend_from_slice(&new_node);
        }

        path_nodes.pop(); // remove trailing comma
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
