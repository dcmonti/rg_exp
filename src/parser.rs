use clap::Parser;

#[derive(Parser, Debug)]
pub struct Args {
    /// Path to the GFA graph file
    #[clap(
        short = 'g',
        long = "graph",
        help = "Input graph (in .gfa format)",
        required = true,
    )]
    pub graph: String,

    /// Path to the reads file
    #[clap(
        short = 'r',
        long = "reads",
        help = "Input reads",
        required = true,
    )]
    pub reads: String,
}
