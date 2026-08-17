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
    about = "Analyze amplicon flank positions and lengths in paired FASTQ files."
)]
struct Args {
    /// R1 FASTQ(.gz)
    #[arg(long)]
    r1: PathBuf,

    /// R2 FASTQ(.gz)
    #[arg(long)]
    r2: PathBuf,

    /// START flank file
    #[arg(long)]
    start_flanks: PathBuf,

    /// END flank file
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


/// Result of classifying canonical R1.
///
/// Expected structure:
///
///     R1:
///
///     START forward -> INSERT -> reverse END
///
/// Possible categories:
///
///     None
///         No START forward flank found.
///
///     R1_START
///         START forward found, but no reverse END found.
///
///     R1_START_revEND
///         START forward and reverse END both found.
///
/// The observed length is:
///
///     - end of reverse END, if reverse END exists
///     - otherwise complete read length
#[derive(Debug)]
struct R1Result {
    category: &'static str,

    start: Option<FlankMatch>,

    reverse_end: Option<FlankMatch>,

    observed_length: usize,
}


/// Result of classifying canonical R2.
///
/// Expected structure:
///
///     R2:
///
///     END forward -> INSERT -> reverse START
///
/// Possible categories:
///
///     None
///         No END forward flank found.
///
///     R2_END
///         END forward found, but no reverse START found.
///
///     R2_END_revSTART
///         END forward and reverse START both found.
///
/// The observed length is:
///
///     - end of reverse START, if reverse START exists
///     - otherwise complete read length
#[derive(Debug)]
struct R2Result {
    category: &'static str,

    end: Option<FlankMatch>,

    reverse_start: Option<FlankMatch>,

    observed_length: usize,
}


/// Histogram of paired R1/R2 lengths.
///
/// Key:
///
///     (pair_category, r1_length, r2_length)
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


/// Counts normalized pair categories.
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


/// Create the pair category.
///
/// R1 is always the canonical R1 after swap correction.
/// R2 is always the canonical R2 after swap correction.
///
/// Examples:
///
///     R1_START + R2_END
///         -> R1_START__R2_END
///
///     R1_START_revEND + R2_END_revSTART
///         -> R1_START_revEND__R2_END_revSTART
///
///     None + R2_END
///         -> None__R2_END
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
    let file =
        File::open(path)
            .with_context(|| {
                format!(
                    "Cannot open flank file {:?}",
                    path
                )
            })?;

    let reader =
        BufReader::new(file);

    let mut flanks =
        Vec::new();

    for (line_no, line) in
        reader.lines().enumerate()
    {
        let line =
            line?;

        let line =
            line.trim();

        if line.is_empty()
            || line.starts_with('#')
        {
            continue;
        }

        let fields:
            Vec<&str> =
            line
                .split_whitespace()
                .collect();

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

        let sequence:
            Vec<u8> =
            sequence
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

        flanks.push(
            Flank {
                name,
                sequence,
            }
        );
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
#[inline]
fn hamming_distance(
    a: &[u8],
    b: &[u8],
    max_distance: usize,
) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let mut distance =
        0usize;

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


/// Find the best match of a set of flanks in a read.
///
/// The caller determines whether the flank sequences are
/// forward or reverse-complemented.
///
/// Ranking:
///
/// 1. Lowest Hamming distance
/// 2. Earliest position
fn find_best_flank(
    sequence: &[u8],
    flanks: &[Flank],
    reverse_complemented: bool,
    max_hamming: usize,
) -> Option<FlankMatch> {
    find_best_flank_from(
        sequence,
        flanks,
        reverse_complemented,
        max_hamming,
        0,
    )
}


/// Find the best flank match at or after `min_position`.
fn find_best_flank_from(
    sequence: &[u8],
    flanks: &[Flank],
    reverse_complemented: bool,
    max_hamming: usize,
    min_position: usize,
) -> Option<FlankMatch> {
    let mut best:
        Option<FlankMatch> =
        None;

    for flank in flanks {
        let oriented_sequence;

        if reverse_complemented {
            oriented_sequence =
                reverse_complement(
                    &flank.sequence
                );
        } else {
            oriented_sequence =
                flank.sequence.clone();
        }

        let k =
            oriented_sequence.len();

        if sequence.len() < k {
            continue;
        }

        if min_position
            > sequence.len() - k
        {
            continue;
        }

        for position in
            min_position
                ..=(sequence.len() - k)
        {
            let candidate =
                &sequence[
                    position
                        ..position + k
                ];

            let distance =
                match hamming_distance(
                    candidate,
                    &oriented_sequence,
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
                            < current
                                .hamming_distance
                            ||
                        (
                            distance
                                == current
                                    .hamming_distance
                            &&
                            position
                                < current.position
                        )
                    }
                };

            if is_better {
                best =
                    Some(
                        FlankMatch {
                            flank_name:
                                flank.name.clone(),

                            position,

                            matched_sequence:
                                candidate.to_vec(),

                            hamming_distance:
                                distance,
                        }
                    );
            }
        }
    }

    best
}


/// Find the best forward START flank.
///
/// This is used specifically for orientation detection and
/// canonical R1 classification.
fn find_forward_start(
    sequence: &[u8],
    start_flanks: &[Flank],
    max_hamming: usize,
) -> Option<FlankMatch> {
    find_best_flank(
        sequence,
        start_flanks,
        false,
        max_hamming,
    )
}


/// Find the best forward END flank.
///
/// This is used specifically for orientation detection and
/// canonical R2 classification.
fn find_forward_end(
    sequence: &[u8],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> Option<FlankMatch> {
    find_best_flank(
        sequence,
        end_flanks,
        false,
        max_hamming,
    )
}


/// Find the reverse-complement END flank.
///
/// This is what we expect at the downstream end of canonical R1.
fn find_reverse_end(
    sequence: &[u8],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> Option<FlankMatch> {
    find_best_flank(
        sequence,
        end_flanks,
        true,
        max_hamming,
    )
}


/// Find the reverse-complement START flank.
///
/// This is what we expect at the downstream end of canonical R2.
fn find_reverse_start(
    sequence: &[u8],
    start_flanks: &[Flank],
    max_hamming: usize,
) -> Option<FlankMatch> {
    find_best_flank(
        sequence,
        start_flanks,
        true,
        max_hamming,
    )
}


/// Determine whether the two input reads are swapped.
///
/// Expected canonical arrangement:
///
///     R1 = START forward
///     R2 = END forward
///
/// Therefore:
///
///     START forward in input R2
///
/// and/or
///
///     END forward in input R1
///
/// indicates that the two reads are swapped.
///
/// We deliberately use the forward biological flanks for this
/// decision, because this mirrors the orientation logic used
/// during demultiplexing.
///
/// Returns:
///
///     true  = input reads must be swapped
///     false = input reads already have canonical orientation
fn reads_are_swapped(
    r1_sequence: &[u8],
    r2_sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> bool {
    let r1_start =
        find_forward_start(
            r1_sequence,
            start_flanks,
            max_hamming,
        )
        .is_some();

    let r2_start =
        find_forward_start(
            r2_sequence,
            start_flanks,
            max_hamming,
        )
        .is_some();

    let r1_end =
        find_forward_end(
            r1_sequence,
            end_flanks,
            max_hamming,
        )
        .is_some();

    /*
     * The canonical orientation is:
     *
     *     R1 = START forward
     *     R2 = END forward
     *
     * Therefore any of these observations means that the
     * reads are reversed:
     *
     *     START forward in R2
     *     END forward in R1
     */
    let swapped =
        r2_start
            || r1_end;

    if swapped {
        eprintln!(
            "Detected swapped read orientation: \
             R1_START={} R2_START={} R1_END={} -> swapping reads",
            r1_start,
            r2_start,
            r1_end,
        );
    }

    swapped
}


/// Classify canonical R1.
///
/// Expected:
///
///     START forward
///
/// followed by:
///
///     insert
///
/// followed optionally by:
///
///     reverse-complement END
///
/// Category:
///
///     None
///     R1_START
///     R1_START_revEND
///
/// Length:
///
///     if reverse END found:
///         end.position + end.length
///
///     otherwise:
///         complete read length
fn classify_r1(
    sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> R1Result {
    let read_length =
        sequence.len();

    let start =
        find_forward_start(
            sequence,
            start_flanks,
            max_hamming,
        );

    match start {
        None => {
            let reverse_end =
                find_reverse_end(
                    sequence,
                    end_flanks,
                    max_hamming,
                );

            R1Result {
                category: "None",
                start: None,
                reverse_end,
                observed_length:
                    read_length,
            }
        }

        Some(start) => {
            let start_pos =
                start.position;

            let start_end =
                start.position
                    + start
                        .matched_sequence
                        .len();

            let reverse_end =
                find_best_flank_from(
                    sequence,
                    end_flanks,
                    true,
                    max_hamming,
                    start_end,
                );

            match reverse_end {
                None => R1Result {
                    category: "R1_START",
                    start: Some(start),
                    reverse_end: None,
                    observed_length:
                        read_length
                            .saturating_sub(
                                start_pos
                            ),
                },

                Some(reverse_end) => R1Result {
                    category:
                        "R1_START_revEND",
                    start: Some(start),
                    reverse_end:
                        Some(reverse_end.clone()),
                    observed_length:
                        reverse_end
                            .position
                            .saturating_sub(
                                start_pos
                            ),
                },
            }
        }
    }
}


/// Classify canonical R2.
///
/// Expected:
///
///     END forward
///
/// followed by:
///
///     insert
///
/// followed optionally by:
///
///     reverse-complement START
///
/// Category:
///
///     None
///     R2_END
///     R2_END_revSTART
///
/// Length:
///
///     if reverse START found:
///         end.position + start.length
///
///     otherwise:
///         complete read length
fn classify_r2(
    sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> R2Result {
    let read_length =
        sequence.len();

    let end =
        find_forward_end(
            sequence,
            end_flanks,
            max_hamming,
        );

    match end {
        None => {
            let reverse_start =
                find_reverse_start(
                    sequence,
                    start_flanks,
                    max_hamming,
                );

            R2Result {
                category: "None",
                end: None,
                reverse_start,
                observed_length:
                    read_length,
            }
        }

        Some(end) => {
            let end_pos =
                end.position;

            let end_end =
                end.position
                    + end
                        .matched_sequence
                        .len();

            let reverse_start =
                find_best_flank_from(
                    sequence,
                    start_flanks,
                    true,
                    max_hamming,
                    end_end,
                );

            match reverse_start {
                None => R2Result {
                    category: "R2_END",
                    end: Some(end),
                    reverse_start: None,
                    observed_length:
                        read_length
                            .saturating_sub(
                                end_pos
                            ),
                },

                Some(reverse_start) => R2Result {
                    category:
                        "R2_END_revSTART",
                    end: Some(end),
                    reverse_start:
                        Some(reverse_start.clone()),
                    observed_length:
                        reverse_start
                            .position
                            .saturating_sub(
                                end_pos
                            ),
                },
            }
        }
    }
}


/// Normalize FASTQ read IDs so R1/R2 can be synchronized.
///
/// Examples:
///
///     read123/1 -> read123
///     read123/2 -> read123
///
///     read123 1:N:0:1 -> read123
fn normalize_read_id(
    id: &[u8],
) -> String {
    let id =
        String::from_utf8_lossy(id);

    id.split_whitespace()
        .next()
        .unwrap_or("")
        .trim_end_matches("/1")
        .trim_end_matches("/2")
        .to_string()
}


/// Write an optional flank match.
///
/// Fields:
///
///     flank name
///     position
///     matched sequence
///     Hamming distance
fn write_optional_match(
    writer:
        &mut csv::Writer<BufWriter<File>>,
    m: Option<&FlankMatch>,
) -> Result<()> {
    match m {
        Some(m) => {
            writer.write_field(
                &m.flank_name
            )?;

            writer.write_field(
                m.position
                    .to_string()
            )?;

            writer.write_field(
                String::from_utf8_lossy(
                    &m.matched_sequence
                )
                .as_bytes()
            )?;

            writer.write_field(
                m.hamming_distance
                    .to_string()
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


/// Compute observed amplicon length from flank positions.
///
/// Defaults:
///
///     start_pos = 0 if start flank is absent
///     end_pos = read_length if end flank is absent
///
/// If an end flank is present but upstream of start, it is ignored
/// and the read end is used as fallback.
fn observed_length_from_positions(
    read_length: usize,
    start: Option<&FlankMatch>,
    end: Option<&FlankMatch>,
) -> usize {
    let start_pos =
        start
            .map(|m| m.position)
            .unwrap_or(0);

    let end_pos =
        match end {
            Some(m)
                if m.position >= start_pos =>
            {
                m.position
            }
            _ => read_length,
        };

    end_pos.saturating_sub(start_pos)
}


/// Process paired FASTQ files.
fn process(
    args: &Args,
    start_flanks: &[Flank],
    end_flanks: &[Flank],
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
     * TSV header.
     */
    writer.write_record([
        "read_id",
        "pair_category",

        "r1_category",

        "r1_start_flank",
        "r1_start_position",
        "r1_start_match",
        "r1_start_hamming",

        "r1_rev_end_flank",
        "r1_rev_end_position",
        "r1_rev_end_match",
        "r1_rev_end_hamming",

        "r1_length",

        "r2_category",

        "r2_end_flank",
        "r2_end_position",
        "r2_end_match",
        "r2_end_hamming",

        "r2_rev_start_flank",
        "r2_rev_start_position",
        "r2_rev_start_match",
        "r2_rev_start_hamming",

        "r2_length",

        "reads_swapped",
    ])?;

    let mut pair_categories =
        PairCategories::new();

    let mut histogram =
        Histogram::new();

    let mut read_count =
        0u64;

    let mut swapped_count =
        0u64;

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
             * Both FASTQ files finished.
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
             * Both records available.
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
                 * Verify that R1 and R2 belong to the
                 * same read pair.
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
                 *
                 * IMPORTANT:
                 *
                 * We do NOT reverse-complement either read.
                 *
                 * If the files are swapped, we simply swap
                 * which sequence is considered canonical R1
                 * and which sequence is considered canonical R2.
                 */
                let input_r1_seq =
                    r1.normalize(false);

                let input_r2_seq =
                    r2.normalize(false);

                /*
                 * ------------------------------------------------
                 * STEP 1: Detect whether the FASTQ pair is swapped.
                 * ------------------------------------------------
                 *
                 * Canonical:
                 *
                 *     R1 = START forward
                 *     R2 = END forward
                 *
                 * Swapped evidence:
                 *
                 *     START forward in input R2
                 *
                 *     OR
                 *
                 *     END forward in input R1
                 */
                let swapped =
                    reads_are_swapped(
                        &input_r1_seq,
                        &input_r2_seq,
                        start_flanks,
                        end_flanks,
                        args.max_hamming,
                    );

                if swapped {
                    swapped_count += 1;
                }

                /*
                 * ------------------------------------------------
                 * STEP 2: Establish canonical R1/R2.
                 * ------------------------------------------------
                 *
                 * No sequence transformation happens here.
                 *
                 * We only choose which physical read belongs
                 * to canonical R1 and canonical R2.
                 */
                let (
                    canonical_r1_seq,
                    canonical_r2_seq,
                ) =
                    if swapped {
                        (
                            &input_r2_seq,
                            &input_r1_seq,
                        )
                    } else {
                        (
                            &input_r1_seq,
                            &input_r2_seq,
                        )
                    };

                /*
                 * ------------------------------------------------
                 * STEP 3: Classify canonical R1.
                 * ------------------------------------------------
                 *
                 * Expected:
                 *
                 *     START forward
                 *
                 * followed optionally by:
                 *
                 *     reverse END
                 */
                let r1_result =
                    classify_r1(
                        canonical_r1_seq,
                        start_flanks,
                        end_flanks,
                        args.max_hamming,
                    );

                /*
                 * ------------------------------------------------
                 * STEP 4: Classify canonical R2.
                 * ------------------------------------------------
                 *
                 * Expected:
                 *
                 *     END forward
                 *
                 * followed optionally by:
                 *
                 *     reverse START
                 */
                let r2_result =
                    classify_r2(
                        canonical_r2_seq,
                        start_flanks,
                        end_flanks,
                        args.max_hamming,
                    );

                let r1_length =
                    observed_length_from_positions(
                        canonical_r1_seq.len(),
                        r1_result.start.as_ref(),
                        r1_result.reverse_end.as_ref(),
                    );

                let r2_length =
                    observed_length_from_positions(
                        canonical_r2_seq.len(),
                        r2_result.end.as_ref(),
                        r2_result.reverse_start.as_ref(),
                    );

                /*
                 * ------------------------------------------------
                 * STEP 5: Pair category.
                 * ------------------------------------------------
                 *
                 * R1 and R2 are now canonical, so the category
                 * has a fixed biological meaning.
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
                 * ------------------------------------------------
                 * STEP 6: Histogram.
                 * ------------------------------------------------
                 *
                 * Store the two lengths belonging to the same
                 * read pair.
                 */
                histogram.add(
                    &pair_category,
                    r1_length,
                    r2_length,
                );

                /*
                 * ------------------------------------------------
                 * STEP 7: Write TSV.
                 * ------------------------------------------------
                 */
                writer.write_field(
                    &r1_id
                )?;

                writer.write_field(
                    &pair_category
                )?;

                /*
                 * -------------------------
                 * Canonical R1
                 * -------------------------
                 */
                writer.write_field(
                    r1_result.category
                )?;

                write_optional_match(
                    &mut writer,
                    r1_result.start.as_ref(),
                )?;

                write_optional_match(
                    &mut writer,
                    r1_result.reverse_end.as_ref(),
                )?;

                writer.write_field(
                    r1_length
                        .to_string()
                )?;

                /*
                 * -------------------------
                 * Canonical R2
                 * -------------------------
                 */
                writer.write_field(
                    r2_result.category
                )?;

                write_optional_match(
                    &mut writer,
                    r2_result.end.as_ref(),
                )?;

                write_optional_match(
                    &mut writer,
                    r2_result.reverse_start.as_ref(),
                )?;

                writer.write_field(
                    r2_length
                        .to_string()
                )?;

                /*
                 * Record whether physical input R1/R2 had
                 * to be swapped.
                 */
                writer.write_field(
                    if swapped {
                        "true"
                    } else {
                        "false"
                    }
                )?;

                /*
                 * Finish TSV row.
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
     * Write pair categories.
     */
    write_pair_categories(
        &args.categories,
        &pair_categories,
    )?;

    eprintln!(
        "Finished: {} read pairs",
        read_count
    );

    eprintln!(
        "Detected swapped orientation in {} read pairs ({:.2}%)",
        swapped_count,
        if read_count > 0 {
            100.0 * swapped_count as f64
                / read_count as f64
        } else {
            0.0
        }
    );

    Ok(())
}


/// Write paired R1/R2 length combinations.
///
/// Output:
///
///     pair_category    r1_length    r2_length    count
///
/// The lengths are canonical:
///
///     r1_length = canonical R1 length
///     r2_length = canonical R2 length
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


/// Write pair-category counts.
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
    ) in entries {
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
        "Loaded {} START flanks and {} END flanks",
        start_flanks.len(),
        end_flanks.len()
    );

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

    /*
     * No separate R1/R2 flank definitions are needed.
     *
     * The biological definitions are:
     *
     *     START
     *     END
     *
     * Orientation handling is explicit:
     *
     *     canonical R1:
     *         START forward
     *         reverse END
     *
     *     canonical R2:
     *         END forward
     *         reverse START
     *
     * Before classification we detect whether the physical
     * FASTQ files have been swapped.
     */
    process(
        &args,
        &start_flanks,
        &end_flanks,
    )
}