## rg_exp
rg_exp is a Rust tool that post-processes minigraph alignments by realigning gaps
between anchors with RecAlign, producing a single GAF alignment per read.

## How it works
1) Run minigraph on the input GFA + reads and parse its chain output.
2) For each read, keep minigraph anchors as exact matches.
3) For gaps between anchors, extract a subgraph and realign the gap with RecAlign.
4) Merge anchors + realigned gaps into a final GAF (printed to stdout).

## Requirements
- Rust toolchain (for building rg_exp and RecAlign)
- minigraph available on PATH (install via conda)

Install minigraph with conda:
```
mamba create -c bioconda -n rg-exps
mamba install -n rg-exps -c bioconda minigraph
```

RecAlign is provided as a git submodule under rg_exp/recalign:
```
git submodule update --init --recursive
```

Build rg_exp and RecAlign:
```
# build rg_exp
cargo build --release

# build RecAlign (submodule)
cd recalign
cargo build --release
```

## Usage
Run with a graph and reads (GFA + FASTA/FASTQ):
```
./target/release/rg_exp -g graph.gfa -r reads.fa -t 8 > alignments.gaf
```

Notes:
- minigraph is invoked internally (hardcoded with -t16).
- If you already have a minigraph chain file, you can skip minigraph:
```
./target/release/rg_exp -g graph.gfa -r reads.fa -c minigraph.chain > alignments.gaf
```
- The thread flag is `-t` / `--threasd` (spelling matches the current CLI).

## Output
- Final GAF alignments are written to stdout.
- Minigraph stdout is mirrored to minigraph_output.txt when minigraph is run.

## Optional environment variables
- RECALIGN_TIMEOUT_SECS: set a timeout (seconds) for each RecAlign gap alignment.
  Use 0 or unset for no timeout.
