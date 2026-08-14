use anyhow::{bail, Context, Result};
use clap::Parser;
use csv::WriterBuilder;
use needletail::{parse_fastx_file, Sequence};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};


/// Command line arguments.
#[derive(Parser, Debug)]
#[command(
    name = "amplicon-qc",
    version,
    about = "Analyze amplicon flank positions and between-flank lengths in paired FASTQ files."
)]
struct Args {
    /// R1 FASTQ(.gz)
    #[arg(long)]
    r1: PathBuf,

    /// R2 FASTQ(.gz)
    #[arg(long)]
    r2: PathBuf,

    /// Start flank file
    #[arg(long)]
    start_flanks: PathBuf,

    /// End flank file
    #[arg(long)]
    end_flanks: PathBuf,

    /// Maximum allowed Hamming distance
    #[arg(long, default_value_t = 1)]
    max_hamming: usize,

    /// Per-read TSV output
    #[arg(long)]
    output: PathBuf,

    /// Aggregated histogram output
    #[arg(long)]
    histogram: PathBuf,

    /// Pair-category count output
    #[arg(long)]
    categories: PathBuf,

    /// Optional maximum number of read pairs to process
    #[arg(long)]
    limit: Option<u64>,
}


/// A flank definition.
#[derive(Debug, Clone)]
struct Flank {
    name: String,
    sequence: Vec<u8>,
}


/// A flank match within a read.
#[derive(Debug, Clone)]
struct FlankMatch {
    flank_name: String,

    /// Zero-based position of the first base of the match.
    position: usize,

    /// Actual sequence found in the read.
    matched_sequence: Vec<u8>,

    /// Hamming distance to the flank definition.
    hamming_distance: usize,
}


/// Classification and length information for one read.
#[derive(Debug)]
struct ReadResult {
    /// One of:
    /// none
    /// start_only
    /// end_only
    /// both
    category: &'static str,

    start: Option<FlankMatch>,
    end: Option<FlankMatch>,

    /// Relevant observed sequence length:
    ///
    /// both:
    ///     sequence between start and end flank
    ///
    /// start_only:
    ///     sequence after start flank until read end
    ///
    /// end_only:
    ///     sequence from read start until end flank
    ///
    /// none:
    ///     complete read
    observed_length: usize,
}


/// Histogram separated by paired-read category.
///
/// Key:
///     (pair_category, length)
///
/// R1 and R2 use separate Histogram objects.
#[derive(Debug)]
struct Histogram {
    lengths: HashMap<(String, usize, usize), u64>,
}

impl Histogram {
    fn new() -> Self {
        Self {
            lengths: HashMap::new(),
        }
    }

    /// Add one complete read-pair.
    ///
    /// The R1 and R2 lengths remain separate.
    ///
    /// Key:
    ///     (pair_category, r1_length, r2_length)
    ///
    /// Example:
    ///     ("end_only__end_only", 23, 22)
    ///
    /// means:
    ///     one read pair where R1 has length 23
    ///     and R2 has length 22.
    fn add(
        &mut self,
        pair_category: &str,
        r1_length: usize,
        r2_length: usize,
    ) {
        *self
            .lengths
            .entry((
                pair_category.to_string(),
                r1_length,
                r2_length,
            ))
            .or_insert(0) += 1;
    }
}


/// Counts the joint R1/R2 categories.
///
/// Example:
///
///     both__both
///     both__start_only
///     start_only__both
///     none__none
#[derive(Debug)]
struct PairCategories {
    counts: HashMap<String, u64>,
}

impl PairCategories {
    fn new() -> Self {
        Self {
            counts: HashMap::new(),
        }
    }

    fn add(
        &mut self,
        r1_category: &'static str,
        r2_category: &'static str,
    ) {
        let category =
            pair_category(
                r1_category,
                r2_category,
            );

        *self
            .counts
            .entry(category)
            .or_insert(0) += 1;
    }
}


/// Generate a combined paired-read category.
fn pair_category(
    r1_category: &'static str,
    r2_category: &'static str,
) -> String {
    format!(
        "{}__{}",
        r1_category,
        r2_category
    )
}


/// Reverse-complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            other => other,
        })
        .collect()
}


/// Read flank definitions.
///
/// Accepted formats:
///
///     ACTGACTG
///
/// or:
///
///     flank_name    ACTGACTG
fn read_flanks(path: &Path) -> Result<Vec<Flank>> {
    let file = File::open(path)
        .with_context(|| {
            format!(
                "Cannot open flank file {:?}",
                path
            )
        })?;

    let reader = BufReader::new(file);

    let mut flanks = Vec::new();

    for (line_no, line) in
        reader.lines().enumerate()
    {
        let line = line?;

        let line = line.trim();

        if line.is_empty()
            || line.starts_with('#')
        {
            continue;
        }

        let fields: Vec<&str> =
            line.split_whitespace().collect();

        let (name, sequence) =
            match fields.len() {
                1 => (
                    format!(
                        "flank_{}",
                        flanks.len()
                    ),
                    fields[0],
                ),

                2 => (
                    fields[0].to_string(),
                    fields[1],
                ),

                _ => bail!(
                    "Invalid flank file {:?}, line {}",
                    path,
                    line_no + 1
                ),
            };

        let sequence: Vec<u8> = sequence
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();

        if sequence.is_empty() {
            bail!(
                "Empty flank sequence in {:?}, line {}",
                path,
                line_no + 1
            );
        }

        flanks.push(Flank {
            name,
            sequence,
        });
    }

    if flanks.is_empty() {
        bail!(
            "No flanks found in {:?}",
            path
        );
    }

    Ok(flanks)
}


/// Return Hamming distance if <= max_distance.
///
/// Returns None immediately once the allowed number
/// of mismatches is exceeded.
#[inline]
fn hamming_distance(
    a: &[u8],
    b: &[u8],
    max_distance: usize,
) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let mut distance = 0usize;

    for i in 0..a.len() {
        if a[i] != b[i] {
            distance += 1;

            if distance > max_distance {
                return None;
            }
        }
    }

    Some(distance)
}


/// Find the best match of any flank in a read.
///
/// Ranking:
///
/// 1. Lowest Hamming distance
/// 2. Earliest position
fn find_flank(
    sequence: &[u8],
    flanks: &[Flank],
    max_hamming: usize,
) -> Option<FlankMatch> {
    let mut best: Option<FlankMatch> = None;

    for flank in flanks {
        let k = flank.sequence.len();

        if sequence.len() < k {
            continue;
        }

        for position in
            0..=(sequence.len() - k)
        {
            let candidate =
                &sequence[position..position + k];

            let distance =
                match hamming_distance(
                    candidate,
                    &flank.sequence,
                    max_hamming,
                ) {
                    Some(d) => d,
                    None => continue,
                };

            let is_better =
                match &best {
                    None => true,

                    Some(current) => {
                        distance
                            < current.hamming_distance
                            ||
                        (
                            distance
                                == current.hamming_distance
                                && position
                                    < current.position
                        )
                    }
                };

            if is_better {
                best = Some(
                    FlankMatch {
                        flank_name:
                            flank.name.clone(),

                        position,

                        matched_sequence:
                            candidate.to_vec(),

                        hamming_distance:
                            distance,
                    },
                );
            }
        }
    }

    best
}


/// Find the best valid START -> END combination.
///
/// Important:
/// We do NOT independently select the best start and end and
/// then check their order. Instead, we explicitly require:
///
///     start_end <= end.position
///
/// The best pair is selected by:
///
/// 1. lowest combined Hamming distance
/// 2. earliest start
/// 3. earliest end
fn find_ordered_flanks(
    sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> Option<(FlankMatch, FlankMatch)> {
    let mut best:
        Option<(FlankMatch, FlankMatch)> = None;

    for start_flank in start_flanks {
        let start_len =
            start_flank.sequence.len();

        if sequence.len() < start_len {
            continue;
        }

        for start_position in
            0..=(sequence.len() - start_len)
        {
            let start_candidate =
                &sequence[
                    start_position
                        ..start_position + start_len
                ];

            let start_distance =
                match hamming_distance(
                    start_candidate,
                    &start_flank.sequence,
                    max_hamming,
                ) {
                    Some(d) => d,
                    None => continue,
                };

            let start_match =
                FlankMatch {
                    flank_name:
                        start_flank.name.clone(),

                    position:
                        start_position,

                    matched_sequence:
                        start_candidate.to_vec(),

                    hamming_distance:
                        start_distance,
                };

            let start_end =
                start_position + start_len;

            /*
             * Search for END only after the START.
             */
            for end_flank in end_flanks {
                let end_len =
                    end_flank.sequence.len();

                if sequence.len() < end_len {
                    continue;
                }

                /*
                 * End has to begin at or after the end
                 * of the start flank.
                 */
                if start_end
                    > sequence.len() - end_len
                {
                    continue;
                }

                for end_position in
                    start_end
                        ..=(sequence.len() - end_len)
                {
                    let end_candidate =
                        &sequence[
                            end_position
                                ..end_position + end_len
                        ];

                    let end_distance =
                        match hamming_distance(
                            end_candidate,
                            &end_flank.sequence,
                            max_hamming,
                        ) {
                            Some(d) => d,
                            None => continue,
                        };

                    let end_match =
                        FlankMatch {
                            flank_name:
                                end_flank.name.clone(),

                            position:
                                end_position,

                            matched_sequence:
                                end_candidate.to_vec(),

                            hamming_distance:
                                end_distance,
                        };

                    let is_better =
                        match &best {
                            None => true,

                            Some((
                                best_start,
                                best_end,
                            )) => {
                                let current_score =
                                    (
                                        start_distance
                                            + end_distance,
                                        start_position,
                                        end_position,
                                    );

                                let best_score =
                                    (
                                        best_start
                                            .hamming_distance
                                            + best_end
                                                .hamming_distance,
                                        best_start.position,
                                        best_end.position,
                                    );

                                current_score
                                    < best_score
                            }
                        };

                    if is_better {
                        best = Some((
                            start_match.clone(),
                            end_match,
                        ));
                    }
                }
            }
        }
    }

    best
}


/// Classify a single read.
///
/// Important behavior:
///
/// - none:
///     length = entire read
///
/// - start_only:
///     length = sequence after start flank
///
/// - end_only:
///     length = sequence before end flank
///
/// - both:
///     only if a valid ordered START -> END pair exists
///
/// If both flank types occur somewhere but no valid ordered
/// combination exists, we fall back to the best start flank
/// and classify as start_only. This keeps the four-category
/// system while ensuring that "both" always means a valid
/// ordered pair.
fn classify_read(
    sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> ReadResult {
    let read_length =
        sequence.len();

    /*
     * First search for a valid ordered pair.
     *
     * This is the only condition under which we call
     * the read "both".
     */
    if let Some((start, end)) =
        find_ordered_flanks(
            sequence,
            start_flanks,
            end_flanks,
            max_hamming,
        )
    {
        let start_end =
            start.position
                + start.matched_sequence.len();

        let observed_length =
            end.position - start_end;

        return ReadResult {
            category: "both",
            start: Some(start),
            end: Some(end),
            observed_length,
        };
    }

    /*
     * No valid START -> END pair.
     *
     * Now independently determine whether either flank
     * exists.
     */
    let start =
        find_flank(
            sequence,
            start_flanks,
            max_hamming,
        );

    let end =
        find_flank(
            sequence,
            end_flanks,
            max_hamming,
        );

    match (start, end) {
        /*
         * No flanks.
         *
         * The entire read is the observed sequence.
         */
        (None, None) => ReadResult {
            category: "none",
            start: None,
            end: None,
            observed_length:
                read_length,
        },

        /*
         * Start only.
         *
         * Everything after the start flank.
         */
        (Some(start), None) => {
            let start_end =
                start.position
                    + start.matched_sequence.len();

            let observed_length =
                read_length
                    .saturating_sub(start_end);

            ReadResult {
                category: "start_only",
                start: Some(start),
                end: None,
                observed_length,
            }
        }

        /*
         * End only.
         *
         * Everything before the end flank.
         */
        (None, Some(end)) => {
            let observed_length =
                end.position;

            ReadResult {
                category: "end_only",
                start: None,
                end: Some(end),
                observed_length,
            }
        }

        /*
         * Both flank types occur, but no ordered pair
         * exists.
         *
         * We cannot call this "both", because "both" is
         * explicitly defined as START -> END.
         *
         * We therefore use the start flank as the
         * anchor and classify as start_only.
         */
        (Some(start), Some(_end)) => {
            let start_end =
                start.position
                    + start.matched_sequence.len();

            let observed_length =
                read_length
                    .saturating_sub(start_end);

            ReadResult {
                category: "start_only",
                start: Some(start),
                end: None,
                observed_length,
            }
        }
    }
}


/// Normalize FASTQ read IDs so R1/R2 can be synchronized.
///
/// For example:
///
/// @read123/1
/// @read123/2
///
/// or:
///
/// @read123 1:N:0:1
///
/// are normalized to the first whitespace-delimited token.
fn normalize_read_id(id: &[u8]) -> String {
    let id =
        String::from_utf8_lossy(id);

    id.split_whitespace()
        .next()
        .unwrap_or("")
        .trim_end_matches("/1")
        .trim_end_matches("/2")
        .to_string()
}


/// Write a flank match to the TSV.
///
/// Writes four fields:
///
/// flank name
/// zero-based position
/// matched sequence
/// Hamming distance
fn write_optional_match(
    writer:
        &mut csv::Writer<BufWriter<File>>,
    m: Option<&FlankMatch>,
) -> Result<()> {
    match m {
        Some(m) => {
            writer.write_field(
                &m.flank_name,
            )?;

            writer.write_field(
                m.position.to_string(),
            )?;

            writer.write_field(
                String::from_utf8_lossy(
                    &m.matched_sequence,
                )
                .as_bytes(),
            )?;

            writer.write_field(
                m.hamming_distance
                    .to_string(),
            )?;
        }

        None => {
            writer.write_field("")?;
            writer.write_field("")?;
            writer.write_field("")?;
            writer.write_field("")?;
        }
    }

    Ok(())
}


/// Process paired FASTQ files.
fn process(
    args: &Args,
    start_flanks_r1: &[Flank],
    end_flanks_r1: &[Flank],
    start_flanks_r2: &[Flank],
    end_flanks_r2: &[Flank],
) -> Result<()> {
    let mut r1_reader =
        parse_fastx_file(&args.r1)
            .with_context(|| {
                format!(
                    "Cannot open {:?}",
                    args.r1
                )
            })?;

    let mut r2_reader =
        parse_fastx_file(&args.r2)
            .with_context(|| {
                format!(
                    "Cannot open {:?}",
                    args.r2
                )
            })?;

    let output_file =
        File::create(&args.output)
            .with_context(|| {
                format!(
                    "Cannot create {:?}",
                    args.output
                )
            })?;

    let mut writer =
        WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(
                BufWriter::new(output_file)
            );

    /*
     * Header.
     */
    writer.write_record([
        "read_id",
        "pair_category",

        "r1_category",
        "r1_start_flank",
        "r1_start_position",
        "r1_start_match",
        "r1_start_hamming",
        "r1_end_flank",
        "r1_end_position",
        "r1_end_match",
        "r1_end_hamming",
        "r1_observed_length",

        "r2_category",
        "r2_start_flank",
        "r2_start_position",
        "r2_start_match",
        "r2_start_hamming",
        "r2_end_flank",
        "r2_end_position",
        "r2_end_match",
        "r2_end_hamming",
        "r2_observed_length",
    ])?;

    /*
     * Joint pair categories.
     */
    let mut pair_categories =
        PairCategories::new();

    /*
     * Separate histograms for R1 and R2.
     *
     * Both are keyed by pair_category + length.
     */
    let mut histogram = Histogram::new();
    let mut read_count = 0u64;

    loop {
        /*
         * Optional processing limit.
         */
        if let Some(limit) =
            args.limit
        {
            if read_count >= limit {
                break;
            }
        }

        let r1_record =
            r1_reader.next();

        let r2_record =
            r2_reader.next();

        match (
            r1_record,
            r2_record,
        ) {
            /*
             * Both files finished normally.
             */
            (None, None) => {
                break;
            }

            /*
             * R1 shorter than R2.
             */
            (None, Some(_)) => {
                bail!(
                    "R1 ended before R2"
                );
            }

            /*
             * R2 shorter than R1.
             */
            (Some(_), None) => {
                bail!(
                    "R2 ended before R1"
                );
            }

            /*
             * Both records exist and parsed correctly.
             */
            (
                Some(Ok(r1)),
                Some(Ok(r2)),
            ) => {
                let r1_id =
                    normalize_read_id(
                        r1.id()
                    );

                let r2_id =
                    normalize_read_id(
                        r2.id()
                    );

                /*
                 * Ensure R1 and R2 are synchronized.
                 */
                if r1_id != r2_id {
                    bail!(
                        "Paired reads out of sync:\n\
                         R1={}\n\
                         R2={}",
                        r1_id,
                        r2_id
                    );
                }

                /*
                 * Normalize sequences to uppercase.
                 */
                let r1_seq =
                    r1.normalize(false);

                let r2_seq =
                    r2.normalize(false);

                /*
                 * Classify independently first.
                 *
                 * They are subsequently combined into
                 * one pair_category.
                 */
                let r1_result =
                    classify_read(
                        &r1_seq,
                        start_flanks_r1,
                        end_flanks_r1,
                        args.max_hamming,
                    );

                let r2_result =
                    classify_read(
                        &r2_seq,
                        start_flanks_r2,
                        end_flanks_r2,
                        args.max_hamming,
                    );

                /*
                 * Joint category.
                 *
                 * Examples:
                 *
                 * both__both
                 * both__start_only
                 * start_only__both
                 * none__none
                 */
                let pair_category =
                    pair_category(
                        r1_result.category,
                        r2_result.category,
                    );

                pair_categories.add(
                    r1_result.category,
                    r2_result.category,
                );

                /*
                 * Add R1 and R2 independently to their
                 * respective histograms, but using the
                 * SAME pair category.
                 *
                 * This is exactly what we want for comparing
                 * R1 and R2 distributions within a paired
                 * category.
                 */
                histogram.add(
                    &pair_category,
                    r1_result.observed_length,
                    r2_result.observed_length,
                );

                /*
                 * Write one TSV row per read pair.
                 */
                writer.write_field(
                    &r1_id,
                )?;

                writer.write_field(
                    &pair_category,
                )?;

                /*
                 * -------------------------
                 * R1
                 * -------------------------
                 */
                writer.write_field(
                    r1_result.category,
                )?;

                write_optional_match(
                    &mut writer,
                    r1_result.start.as_ref(),
                )?;

                write_optional_match(
                    &mut writer,
                    r1_result.end.as_ref(),
                )?;

                writer.write_field(
                    r1_result
                        .observed_length
                        .to_string(),
                )?;

                /*
                 * -------------------------
                 * R2
                 * -------------------------
                 */
                writer.write_field(
                    r2_result.category,
                )?;

                write_optional_match(
                    &mut writer,
                    r2_result.start.as_ref(),
                )?;

                write_optional_match(
                    &mut writer,
                    r2_result.end.as_ref(),
                )?;

                writer.write_field(
                    r2_result
                        .observed_length
                        .to_string(),
                )?;

                /*
                 * Finish TSV row.
                 *
                 * IMPORTANT:
                 * There is intentionally only ONE
                 * write_record() here.
                 */
                writer.write_record(
                    None::<&[u8]>
                )?;

                read_count += 1;

                if read_count % 1_000_000 == 0 {
                    eprintln!(
                        "Processed {:>12} read pairs",
                        read_count
                    );
                }
            }

            /*
             * FASTQ parsing error in R1.
             */
            (Some(Err(e)), _) => {
                return Err(e.into());
            }

            /*
             * FASTQ parsing error in R2.
             */
            (_, Some(Err(e))) => {
                return Err(e.into());
            }
        }
    }

    writer.flush()?;

    /*
     * Write histogram.
     */
    write_histogram(
        &args.histogram,
        &histogram,
    )?;

    /*
     * Write joint category counts.
     */
    write_pair_categories(
        &args.categories,
        &pair_categories,
    )?;

    eprintln!(
        "Finished: {} read pairs",
        read_count
    );

    Ok(())
}


/// Write paired R1/R2 length combinations.
///
/// Output:
///
/// pair_category    r1_length    r2_length    count
///
/// Each row represents one observed R1/R2 length combination
/// within a paired-read category.
///
/// Example:
///
/// end_only__end_only    23    22    270
/// end_only__end_only    23    23     15
/// end_only__end_only    24    22     32
///
/// This does NOT add R1 and R2 lengths.
///
/// Instead, it preserves the fact that the two lengths
/// belonged to the same read pair.
fn write_histogram(
    path: &Path,
    histogram: &Histogram,
) -> Result<()> {
    let file =
        File::create(path)?;

    let mut writer =
        BufWriter::new(file);

    writeln!(
        writer,
        "pair_category\tr1_length\tr2_length\tcount"
    )?;

    let mut entries:
        Vec<_> =
        histogram
            .lengths
            .iter()
            .collect();

    entries.sort_unstable_by(
        |(
            (cat_a, r1_a, r2_a),
            _,
        ),
         (
            (cat_b, r1_b, r2_b),
            _,
         )| {
            cat_a
                .cmp(cat_b)
                .then_with(|| r1_a.cmp(r1_b))
                .then_with(|| r2_a.cmp(r2_b))
        },
    );

    for (
        (pair_category, r1_length, r2_length),
        count,
    ) in entries {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            pair_category,
            r1_length,
            r2_length,
            count,
        )?;
    }

    writer.flush()?;

    Ok(())
}

/// Write counts for the joint R1/R2 categories.
fn write_pair_categories(
    path: &Path,
    categories: &PairCategories,
) -> Result<()> {
    let file =
        File::create(path)?;

    let mut writer =
        BufWriter::new(file);

    writeln!(
        writer,
        "pair_category\tcount"
    )?;

    let mut entries:
        Vec<_> =
        categories
            .counts
            .iter()
            .collect();

    entries.sort_unstable_by(
        |(category_a, _),
         (category_b, _)| {
            category_a.cmp(category_b)
        },
    );

    for (
        category,
        count,
    ) in entries
    {
        writeln!(
            writer,
            "{}\t{}",
            category,
            count,
        )?;
    }

    writer.flush()?;

    Ok(())
}


/// Program entry point.
fn main() -> Result<()> {
    let args =
        Args::parse();

    eprintln!(
        "Reading flank definitions..."
    );

    let start_flanks =
        read_flanks(
            &args.start_flanks
        )?;

    let end_flanks =
        read_flanks(
            &args.end_flanks
        )?;

    eprintln!(
        "Loaded {} start flanks and {} end flanks",
        start_flanks.len(),
        end_flanks.len()
    );

    /*
     * R1 orientation:
     *
     * START -> AMPLICON -> END
     *
     *
     * R2 orientation:
     *
     * revcomp(END)
     *      ->
     * revcomp(AMPLICON)
     *      ->
     * revcomp(START)
     *
     * Therefore:
     *
     * R2 START = reverse-complement of R1 END
     * R2 END   = reverse-complement of R1 START
     */

    let start_flanks_r2:
        Vec<Flank> =
        end_flanks
            .iter()
            .map(|flank| {
                Flank {
                    name:
                        flank.name.clone(),

                    sequence:
                        reverse_complement(
                            &flank.sequence
                        ),
                }
            })
            .collect();

    let end_flanks_r2:
        Vec<Flank> =
        start_flanks
            .iter()
            .map(|flank| {
                Flank {
                    name:
                        flank.name.clone(),

                    sequence:
                        reverse_complement(
                            &flank.sequence
                        ),
                }
            })
            .collect();

    eprintln!(
        "Maximum Hamming distance: {}",
        args.max_hamming
    );

    if let Some(limit) =
        args.limit
    {
        eprintln!(
            "Read-pair limit: {}",
            limit
        );
    }

    process(
        &args,
        &start_flanks,
        &end_flanks,
        &start_flanks_r2,
        &end_flanks_r2,
    )
}