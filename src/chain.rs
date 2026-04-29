use std::fmt::Display;

use gfa::gfa::GFA;
use crate::graph_index;

// Walks forward through `nodes` starting at `start_from` until accumulated length
// (starting from `already_accumulated`) reaches `target`.
// Returns (last_idx, trim): last_idx is the last included node index; trim is Some(k) when the
// boundary node should keep only its first k bytes.
fn walk_forward(
    nodes: &[(Vec<u8>, bool)],
    index: &graph_index::PathIndex,
    start_from: usize,
    already_accumulated: usize,
    target: usize,
) -> (usize, Option<usize>) {
    if already_accumulated >= target || start_from >= nodes.len() {
        return (start_from.saturating_sub(1), None);
    }
    let mut curr = already_accumulated;
    let mut en = start_from;
    let mut trim = None;
    while curr < target && en < nodes.len() {
        let nlen = index.get(&nodes[en].0).unwrap().segment.len();
        if curr + nlen >= target {
            trim = Some(target - curr);
            curr = target;
        } else {
            curr += nlen;
        }
        en += 1;
    }
    (en.saturating_sub(1), trim)
}

// Walks backward through `nodes` starting just below `anchor_idx` until accumulated length
// (starting from `already_accumulated`) reaches `target`.
// Returns (first_idx, trim): first_idx is the lowest included index; trim is Some(k) when the
// boundary node should keep only its last k bytes.
fn walk_backward(
    nodes: &[(Vec<u8>, bool)],
    index: &graph_index::PathIndex,
    anchor_idx: usize,
    already_accumulated: usize,
    target: usize,
) -> (usize, Option<usize>) {
    if already_accumulated >= target {
        return (anchor_idx, None);
    }
    let mut curr = already_accumulated;
    let mut st = anchor_idx;
    let mut trim = None;
    while curr < target && st > 0 {
        st -= 1;
        let nlen = index.get(&nodes[st].0).unwrap().segment.len();
        if curr + nlen >= target {
            trim = Some(target - curr);
            curr = target;
        } else {
            curr += nlen;
        }
    }
    (st, trim)
}

pub fn build_chain(alignment: &str) -> Chain {
    // Chain from minigraph alignment output parsing logic goes here
    // Expected input:
    // * sName sLen nMz div sStart sEnd qStart qEnd\n
    let mut chain = Chain::new(vec![]);
    alignment.lines().for_each(|line| {
        let fields = line.split('\t').collect::<Vec<&str>>();
        if fields.len() == 9 {
            let s_name = &fields[1][1..];
            let s_len: usize = fields[2].parse().unwrap();
            let mut s_start: usize = fields[5].parse().unwrap();
            let mut s_end: usize = fields[6].parse().unwrap();
            let q_start: usize = fields[7].parse().unwrap();
            let q_end: usize = fields[8].parse().unwrap();
            let orientation = fields[1].chars().next().unwrap() == '>';

            if !orientation {
                let tmp = s_start;
                s_start = s_len - s_end;
                s_end = s_len - tmp;
            }

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
    //let paths_to_keep = path_index.both_paths(from_node, to_node);
    let mut paths_to_keep = path_index.common_paths(from_node, to_node);
    if paths_to_keep.len() == 0 {
        eprintln!(
            "Warning: no common paths found between from_node {} and to_node {}. Falling back to both_paths.",
            String::from_utf8_lossy(from_node),
            String::from_utf8_lossy(to_node)
        );
        paths_to_keep = path_index.both_paths(from_node, to_node);
    }
    let mut subgraph: GFA<Vec<u8>, ()> = GFA::new();
    let mut segments_set = std::collections::HashSet::new();

    for (path_id, relation) in paths_to_keep {
        let path_name = &path_index.paths_to_ids[&path_id];
        let path_nodes_list = &path_index.paths[&path_id].nodes;

        // Compute (st_idx, en_idx, boundary_trim) for the subgraph range.
        // boundary_trim = Some((is_low_boundary, keep_len)):
        //   is_low_boundary=true  -> st_idx is the boundary node; keep its last keep_len bytes
        //   is_low_boundary=false -> en_idx is the boundary node; keep its first keep_len bytes
        // Anchor effective lengths account for the partial inclusion of from_node/to_node.
        let (st_idx, en_idx, boundary_trim): (usize, usize, Option<(bool, usize)>) = match relation {
            graph_index::PathRelation::Both => {
                let from_idx = path_index.get(from_node).unwrap().get_path_fl(&path_id).unwrap().0;
                let to_idx   = path_index.get(to_node).unwrap().get_path_fl(&path_id).unwrap().1;
                let (s, e) = if from_idx <= to_idx { (from_idx, to_idx) } else { (to_idx, from_idx) };
                (s, e, None)
            },
            graph_index::PathRelation::First => {
                let anchor_idx = path_index.get(from_node).unwrap().get_path_fl(&path_id).unwrap().0;
                let from_full_len = path_index.get(from_node).unwrap().segment.len();
                if reversed {
                    // from_node kept as [..position.1]; effective contribution = position.1
                    let effective = from_anchor.graph_pos.position.1;
                    let (st, trim) = walk_backward(path_nodes_list, path_index, anchor_idx, effective, substr_len);
                    (st, anchor_idx, trim.map(|k| (true, k)))
                } else {
                    // from_node kept as [position.0..]; effective contribution = full_len - position.0
                    let effective = from_full_len.saturating_sub(from_anchor.graph_pos.position.0);
                    let (en, trim) = walk_forward(path_nodes_list, path_index, anchor_idx + 1, effective, substr_len);
                    (anchor_idx, en, trim.map(|k| (false, k)))
                }
            },
            graph_index::PathRelation::Second => {
                let anchor_idx = path_index.get(to_node).unwrap().get_path_fl(&path_id).unwrap().1;
                let to_full_len = path_index.get(to_node).unwrap().segment.len();
                if reversed {
                    // to_node kept as [position.0..]; effective contribution = full_len - position.0
                    let effective = to_full_len.saturating_sub(to_anchor.graph_pos.position.0);
                    let (en, trim) = walk_forward(path_nodes_list, path_index, anchor_idx + 1, effective, substr_len);
                    (anchor_idx, en, trim.map(|k| (false, k)))
                } else {
                    // to_node kept as [..position.1]; effective contribution = position.1
                    let effective = to_anchor.graph_pos.position.1;
                    let (st, trim) = walk_backward(path_nodes_list, path_index, anchor_idx, effective, substr_len);
                    (st, anchor_idx, trim.map(|k| (true, k)))
                }
            },
        };

        let st_idx = st_idx.min(path_nodes_list.len().saturating_sub(1));
        let en_idx = en_idx.min(path_nodes_list.len().saturating_sub(1));
        let (st_idx, en_idx) = if st_idx <= en_idx { (st_idx, en_idx) } else { (en_idx, st_idx) };
        let range_len = en_idx - st_idx;

        let mut path_nodes = Vec::new();

        // Always iterate in path order. For reversed reads, from_node sits at the high
        // end and to_node at the low end of the range; the RC'd read (sent to recalign)
        // aligns left-to-right against the forward-orientation path segments.
        for (rel_i, (node_id, path_orientation)) in path_nodes_list[st_idx..=en_idx].iter().enumerate() {
            let mut seq_slice: &[u8] = path_index.get(node_id).unwrap().segment.as_slice();

            if node_id == from_node {
                if reversed {
                    // gap starts before anchor start in forward coords: keep [..position.0]
                    seq_slice = &seq_slice[..from_anchor.graph_pos.position.0];
                } else {
                    // gap starts after anchor end in forward coords: keep [position.1..]
                    seq_slice = &seq_slice[from_anchor.graph_pos.position.1..];
                }
            }
            if node_id == to_node {
                if reversed {
                    // gap ends after anchor end in forward coords: keep [position.1..]
                    seq_slice = &seq_slice[to_anchor.graph_pos.position.1..];
                } else {
                    // gap ends before anchor start in forward coords: keep [..position.0]
                    seq_slice = &seq_slice[..to_anchor.graph_pos.position.0];
                }
            }

            if seq_slice.is_empty() {
                continue;
            }
            // Trim the boundary node (only applies to non-anchor nodes in First/Second cases).
            if node_id != from_node && node_id != to_node {
                if let Some((is_low_boundary, keep_len)) = boundary_trim {
                    let is_boundary = (is_low_boundary && rel_i == 0)
                        || (!is_low_boundary && rel_i == range_len);
                    if is_boundary {
                        let len = seq_slice.len();
                        if is_low_boundary {
                            // backward walk: boundary is first in path order; keep its tail
                            seq_slice = &seq_slice[len.saturating_sub(keep_len)..];
                        } else {
                            // forward walk: boundary is last in path order; keep its head
                            seq_slice = &seq_slice[..keep_len.min(len)];
                        }
                    }
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
