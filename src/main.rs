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

    /// Biological START flank file
    #[arg(long)]
    start_flanks: PathBuf,

    /// Biological END flank file
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


/// A biological flank definition.
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

    /// Hamming distance to the biological flank sequence
    /// in the searched orientation.
    hamming_distance: usize,
}


/// Classification of one read based on biological flank identity.
///
/// The classification is independent of whether the read is
/// physically R1 or R2.
///
/// A biological START flank is searched in both orientations:
///
///     START
///     reverse-complement(START)
///
/// A biological END flank is also searched in both orientations:
///
///     END
///     reverse-complement(END)
///
/// A complete read can therefore appear in either physical
/// orientation:
///
///     START -> insert -> END
///
/// or:
///
///     END -> insert -> START
///
/// Both are classified as:
///
///     both
///
/// This is essential because R1 and R2 can be swapped and because
/// a read itself may represent either strand.
#[derive(Debug)]
struct ReadResult {
    category: &'static str,

    /// Biological START flank, if found.
    start: Option<FlankMatch>,

    /// Biological END flank, if found.
    end: Option<FlankMatch>,

    /// Orientation in which the biological START flank
    /// was actually found in the read.
    ///
    /// "forward"
    /// "reverse_complement"
    start_orientation: Option<&'static str>,

    /// Orientation in which the biological END flank
    /// was actually found in the read.
    end_orientation: Option<&'static str>,

    /// Physical position of the biological START flank.
    start_position: Option<usize>,

    /// Physical position of the biological END flank.
    end_position: Option<usize>,

    /// Length of sequence between the two biological flanks
    /// when both are found.
    ///
    /// For single-flank reads:
    ///
    /// start_only:
    ///     sequence after START if START is present
    ///
    /// end_only:
    ///     sequence before END if END is present
    ///
    /// none:
    ///     complete read length
    observed_length: usize,
}


/// Histogram of paired R1/R2 lengths.
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
        r1: &ReadResult,
        r2: &ReadResult,
    ) {
        let category =
            normalized_pair_category(
                r1.category,
                r2.category,
            );

        *self
            .counts
            .entry(category)
            .or_insert(0) += 1;
    }
}


/// Normalize the pair category so that swapping R1 and R2
/// does not create a different biological category.
///
/// Examples:
///
///     start_only + end_only
///         -> start_only__end_only
///
///     end_only + start_only
///         -> start_only__end_only
///
///     both + start_only
///         -> both__start_only
///
///     start_only + both
///         -> both__start_only
///
///     none + end_only
///         -> end_only__none
///
///     none + none
///         -> none__none
fn normalized_pair_category(
    category_a: &'static str,
    category_b: &'static str,
) -> String {
    fn rank(category: &str) -> usize {
        match category {
            "both" => 0,
            "start_only" => 1,
            "end_only" => 2,
            "none" => 3,
            _ => 99,
        }
    }

    if rank(category_a) <= rank(category_b) {
        format!(
            "{}__{}",
            category_a,
            category_b
        )
    } else {
        format!(
            "{}__{}",
            category_b,
            category_a
        )
    }
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
/// Returns None once the allowed number of mismatches
/// is exceeded.
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


/// Find the best match of one biological flank in either
/// orientation.
///
/// Returns:
///
///     FlankMatch
///     orientation
///
/// Ranking:
///
/// 1. Lowest Hamming distance
/// 2. Earliest position
fn find_best_oriented_flank(
    sequence: &[u8],
    flanks: &[Flank],
    max_hamming: usize,
) -> Option<(FlankMatch, &'static str)> {
    let mut best:
        Option<(FlankMatch, &'static str)> =
        None;

    for flank in flanks {
        let reverse =
            reverse_complement(
                &flank.sequence
            );

        let orientations: [
            (&[u8], &'static str);
            2
        ] = [
            (
                &flank.sequence,
                "forward",
            ),
            (
                &reverse,
                "reverse_complement",
            ),
        ];

        for (
            oriented_sequence,
            orientation,
        ) in orientations
        {
            let k =
                oriented_sequence.len();

            if sequence.len() < k {
                continue;
            }

            for position in
                0..=(sequence.len() - k)
            {
                let candidate =
                    &sequence[
                        position
                            ..position + k
                    ];

                let distance =
                    match hamming_distance(
                        candidate,
                        oriented_sequence,
                        max_hamming,
                    ) {
                        Some(d) => d,
                        None => continue,
                    };

                let match_record =
                    FlankMatch {
                        flank_name:
                            flank.name.clone(),

                        position,

                        matched_sequence:
                            candidate.to_vec(),

                        hamming_distance:
                            distance,
                    };

                let is_better =
                    match &best {
                        None => true,

                        Some((
                            current,
                            _,
                        )) => {
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
                                    < current
                                        .position
                            )
                        }
                    };

                if is_better {
                    best = Some((
                        match_record,
                        orientation,
                    ));
                }
            }
        }
    }

    best
}


/// Find a valid COMPLETE biological flank pair.
///
/// IMPORTANT:
///
/// There are TWO valid physical orientations:
///
///     START -> END
///
/// and:
///
///     END -> START
///
/// Both represent a complete biological amplicon and are
/// therefore classified as "both".
///
/// The returned fields are always biological:
///
///     start = biological START flank
///     end   = biological END flank
///
/// regardless of their physical order in the read.
///
/// Ranking:
///
/// 1. Lowest combined Hamming distance
/// 2. Smallest distance between the flanks
/// 3. Earliest first flank
/// 4. Earliest second flank
fn find_ordered_biological_flanks(
    sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> Option<(
    FlankMatch,
    &'static str,
    FlankMatch,
    &'static str,
)> {
    let mut best:
        Option<(
            FlankMatch,
            &'static str,
            FlankMatch,
            &'static str,
        )> = None;

    /*
     * --------------------------------------------------------
     * Case 1:
     *
     * Biological START -> biological END
     * --------------------------------------------------------
     */

    for start_flank in start_flanks {
        let start_reverse =
            reverse_complement(
                &start_flank.sequence
            );

        let start_orientations: [
            (&[u8], &'static str);
            2
        ] = [
            (
                &start_flank.sequence,
                "forward",
            ),
            (
                &start_reverse,
                "reverse_complement",
            ),
        ];

        for (
            start_sequence,
            start_orientation,
        ) in start_orientations
        {
            let start_len =
                start_sequence.len();

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
                        start_sequence,
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
                    start_position
                        + start_len;

                for end_flank in end_flanks {
                    let end_reverse =
                        reverse_complement(
                            &end_flank.sequence
                        );

                    let end_orientations: [
                        (&[u8], &'static str);
                        2
                    ] = [
                        (
                            &end_flank.sequence,
                            "forward",
                        ),
                        (
                            &end_reverse,
                            "reverse_complement",
                        ),
                    ];

                    for (
                        end_sequence,
                        end_orientation,
                    ) in end_orientations
                    {
                        let end_len =
                            end_sequence.len();

                        if sequence.len()
                            < end_len
                        {
                            continue;
                        }

                        if start_end
                            > sequence.len()
                                - end_len
                        {
                            continue;
                        }

                        for end_position in
                            start_end
                                ..=(sequence.len()
                                    - end_len)
                        {
                            let end_candidate =
                                &sequence[
                                    end_position
                                        ..end_position
                                            + end_len
                                ];

                            let end_distance =
                                match hamming_distance(
                                    end_candidate,
                                    end_sequence,
                                    max_hamming,
                                ) {
                                    Some(d) => d,
                                    None => continue,
                                };

                            let end_match =
                                FlankMatch {
                                    flank_name:
                                        end_flank
                                            .name
                                            .clone(),

                                    position:
                                        end_position,

                                    matched_sequence:
                                        end_candidate
                                            .to_vec(),

                                    hamming_distance:
                                        end_distance,
                                };

                            let gap =
                                end_position
                                    - start_end;

                            let current_score = (
                                start_distance
                                    + end_distance,
                                gap,
                                start_position,
                                end_position,
                            );

                            let is_better =
                                match &best {
                                    None => true,

                                    Some((
                                        best_start,
                                        _,
                                        best_end,
                                        _,
                                    )) => {
                                        let best_start_end =
                                            best_start.position
                                                + best_start
                                                    .matched_sequence
                                                    .len();

                                        let best_gap =
                                            best_end.position
                                                .saturating_sub(
                                                    best_start_end
                                                );

                                        let best_score = (
                                            best_start
                                                .hamming_distance
                                                + best_end
                                                    .hamming_distance,
                                            best_gap,
                                            best_start
                                                .position,
                                            best_end
                                                .position,
                                        );

                                        current_score
                                            < best_score
                                    }
                                };

                            if is_better {
                                best = Some((
                                    start_match.clone(),
                                    start_orientation,
                                    end_match,
                                    end_orientation,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     * --------------------------------------------------------
     * Case 2:
     *
     * Biological END -> biological START
     *
     * This is the reverse physical orientation of the same
     * biological amplicon.
     *
     * We return the fields in biological order:
     *
     *     start = START match
     *     end   = END match
     *
     * even though END occurs physically first in the read.
     * --------------------------------------------------------
     */

    for end_flank in end_flanks {
        let end_reverse =
            reverse_complement(
                &end_flank.sequence
            );

        let end_orientations: [
            (&[u8], &'static str);
            2
        ] = [
            (
                &end_flank.sequence,
                "forward",
            ),
            (
                &end_reverse,
                "reverse_complement",
            ),
        ];

        for (
            end_sequence,
            end_orientation,
        ) in end_orientations
        {
            let end_len =
                end_sequence.len();

            if sequence.len() < end_len {
                continue;
            }

            for end_position in
                0..=(sequence.len() - end_len)
            {
                let end_candidate =
                    &sequence[
                        end_position
                            ..end_position + end_len
                    ];

                let end_distance =
                    match hamming_distance(
                        end_candidate,
                        end_sequence,
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

                let end_end =
                    end_position
                        + end_len;

                for start_flank in start_flanks {
                    let start_reverse =
                        reverse_complement(
                            &start_flank.sequence
                        );

                    let start_orientations: [
                        (&[u8], &'static str);
                        2
                    ] = [
                        (
                            &start_flank.sequence,
                            "forward",
                        ),
                        (
                            &start_reverse,
                            "reverse_complement",
                        ),
                    ];

                    for (
                        start_sequence,
                        start_orientation,
                    ) in start_orientations
                    {
                        let start_len =
                            start_sequence.len();

                        if sequence.len()
                            < start_len
                        {
                            continue;
                        }

                        if end_end
                            > sequence.len()
                                - start_len
                        {
                            continue;
                        }

                        for start_position in
                            end_end
                                ..=(sequence.len()
                                    - start_len)
                        {
                            let start_candidate =
                                &sequence[
                                    start_position
                                        ..start_position
                                            + start_len
                                ];

                            let start_distance =
                                match hamming_distance(
                                    start_candidate,
                                    start_sequence,
                                    max_hamming,
                                ) {
                                    Some(d) => d,
                                    None => continue,
                                };

                            let start_match =
                                FlankMatch {
                                    flank_name:
                                        start_flank
                                            .name
                                            .clone(),

                                    position:
                                        start_position,

                                    matched_sequence:
                                        start_candidate
                                            .to_vec(),

                                    hamming_distance:
                                        start_distance,
                                };

                            let gap =
                                start_position
                                    - end_end;

                            let current_score = (
                                end_distance
                                    + start_distance,
                                gap,
                                end_position,
                                start_position,
                            );

                            let is_better =
                                match &best {
                                    None => true,

                                    Some((
                                        best_start,
                                        _,
                                        best_end,
                                        _,
                                    )) => {
                                        let best_first_position =
                                            best_start.position
                                                .min(
                                                    best_end.position
                                                );

                                        let best_last_position =
                                            best_start.position
                                                .max(
                                                    best_end.position
                                                );

                                        let best_start_end =
                                            best_end.position
                                                + best_end
                                                    .matched_sequence
                                                    .len();

                                        let best_gap =
                                            if best_end.position
                                                < best_start.position
                                            {
                                                best_start.position
                                                    .saturating_sub(
                                                        best_start_end
                                                    )
                                            } else {
                                                best_start
                                                    .position
                                                    .saturating_sub(
                                                        best_end.position
                                                            + best_end
                                                                .matched_sequence
                                                                .len()
                                                    )
                                            };

                                        let best_score = (
                                            best_start
                                                .hamming_distance
                                                + best_end
                                                    .hamming_distance,
                                            best_gap,
                                            best_first_position,
                                            best_last_position,
                                        );

                                        current_score
                                            < best_score
                                    }
                                };

                            if is_better {
                                /*
                                 * IMPORTANT:
                                 *
                                 * Return biological START first,
                                 * biological END second.
                                 *
                                 * Even though END physically occurs
                                 * before START in the read.
                                 */
                                best = Some((
                                    start_match,
                                    start_orientation,
                                    end_match,
                                    end_orientation,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    best
}


/// Classify one read using biological flank identity.
///
/// Rules:
///
///     START -> END = both
///     END   -> START = both
///
/// Thus both physical orientations of an amplicon are treated
/// identically.
///
/// The biological START and END matches are always returned in
/// their biological roles, independent of physical order.
fn classify_read(
    sequence: &[u8],
    start_flanks: &[Flank],
    end_flanks: &[Flank],
    max_hamming: usize,
) -> ReadResult {
    let read_length =
        sequence.len();

    /*
     * --------------------------------------------------------
     * First: look for a complete biological pair.
     *
     * This now recognizes BOTH:
     *
     *     START -> END
     *
     * and:
     *
     *     END -> START
     * --------------------------------------------------------
     */

    if let Some((
        start,
        start_orientation,
        end,
        end_orientation,
    )) =
        find_ordered_biological_flanks(
            sequence,
            start_flanks,
            end_flanks,
            max_hamming,
        )
    {
        let start_end =
            start.position
                + start.matched_sequence.len();

        let end_end =
            end.position
                + end.matched_sequence.len();

        let observed_length =
            if start.position < end.position {
                /*
                 * Physical order:
                 *
                 * START -> END
                 */
                end.position
                    .saturating_sub(
                        start_end
                    )
            } else {
                /*
                 * Physical order:
                 *
                 * END -> START
                 */
                start.position
                    .saturating_sub(
                        end_end
                    )
            };

        return ReadResult {
            category: "both",

            start_position:
                Some(start.position),

            end_position:
                Some(end.position),

            start: Some(start),

            end: Some(end),

            start_orientation:
                Some(start_orientation),

            end_orientation:
                Some(end_orientation),

            observed_length,
        };
    }

    /*
     * --------------------------------------------------------
     * No complete biological pair.
     *
     * Search independently for START and END.
     * --------------------------------------------------------
     */

    let start =
        find_best_oriented_flank(
            sequence,
            start_flanks,
            max_hamming,
        );

    let end =
        find_best_oriented_flank(
            sequence,
            end_flanks,
            max_hamming,
        );

    match (start, end) {
        /*
         * ----------------------------------------------------
         * NONE
         * ----------------------------------------------------
         */
        (None, None) => ReadResult {
            category: "none",

            start: None,
            end: None,

            start_orientation: None,
            end_orientation: None,

            start_position: None,
            end_position: None,

            observed_length:
                read_length,
        },

        /*
         * ----------------------------------------------------
         * START ONLY
         * ----------------------------------------------------
         */
        (
            Some((
                start,
                start_orientation,
            )),
            None,
        ) => {
            let start_end =
                start.position
                    + start.matched_sequence.len();

            let observed_length =
                read_length
                    .saturating_sub(
                        start_end
                    );

            ReadResult {
                category: "start_only",

                start_position:
                    Some(start.position),

                end_position: None,

                start: Some(start),
                end: None,

                start_orientation:
                    Some(start_orientation),

                end_orientation: None,

                observed_length,
            }
        }

        /*
         * ----------------------------------------------------
         * END ONLY
         * ----------------------------------------------------
         */
        (
            None,
            Some((
                end,
                end_orientation,
            )),
        ) => {
            let observed_length =
                end.position;

            ReadResult {
                category: "end_only",

                start_position: None,

                end_position:
                    Some(end.position),

                start: None,
                end: Some(end),

                start_orientation: None,

                end_orientation:
                    Some(end_orientation),

                observed_length,
            }
        }

        /*
         * ----------------------------------------------------
         * BOTH FLANK TYPES FOUND, BUT NO VALID ORDERED PAIR
         * ----------------------------------------------------
         *
         * This is an ambiguous / conflicting situation.
         *
         * We do NOT silently discard the END flank anymore.
         *
         * Instead we choose the biological flank pair with
         * the best independent matches and classify according
         * to which flank is the better anchor.
         *
         * In practice this should be rare if the flank
         * definitions are sufficiently specific.
         * ----------------------------------------------------
         */
        (
            Some((
                start,
                start_orientation,
            )),
            Some((
                end,
                end_orientation,
            )),
        ) => {
            let start_quality =
                start.hamming_distance;

            let end_quality =
                end.hamming_distance;

            if start_quality <= end_quality {
                let start_end =
                    start.position
                        + start.matched_sequence.len();

                let observed_length =
                    read_length
                        .saturating_sub(
                            start_end
                        );

                ReadResult {
                    category: "start_only",

                    start_position:
                        Some(start.position),

                    end_position: None,

                    start: Some(start),
                    end: None,

                    start_orientation:
                        Some(start_orientation),

                    end_orientation: None,

                    observed_length,
                }
            } else {
                let observed_length =
                    end.position;

                ReadResult {
                    category: "end_only",

                    start_position: None,

                    end_position:
                        Some(end.position),

                    start: None,
                    end: Some(end),

                    start_orientation: None,

                    end_orientation:
                        Some(end_orientation),

                    observed_length,
                }
            }
        }
    }
}


/// Normalize FASTQ read IDs so R1/R2 can be synchronized.
///
/// Examples:
///
/// @read123/1
/// @read123/2
///
/// become:
///
/// read123
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


/// Write an optional flank match to the TSV.
fn write_optional_match(
    writer:
        &mut csv::Writer<BufWriter<File>>,
    m: Option<&FlankMatch>,
    orientation: Option<&'static str>,
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

            writer.write_field(
                orientation.unwrap_or(""),
            )?;
        }

        None => {
            writer.write_field("")?;
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
        "r1_start_orientation",
        "r1_end_flank",
        "r1_end_position",
        "r1_end_match",
        "r1_end_hamming",
        "r1_end_orientation",
        "r1_observed_length",

        "r2_category",
        "r2_start_flank",
        "r2_start_position",
        "r2_start_match",
        "r2_start_hamming",
        "r2_start_orientation",
        "r2_end_flank",
        "r2_end_position",
        "r2_end_match",
        "r2_end_hamming",
        "r2_end_orientation",
        "r2_observed_length",
    ])?;

    let mut pair_categories =
        PairCategories::new();

    let mut histogram =
        Histogram::new();

    let mut read_count =
        0u64;

    loop {
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
             * Both files ended normally.
             */
            (None, None) => {
                break;
            }

            /*
             * R1 ended too early.
             */
            (None, Some(_)) => {
                bail!(
                    "R1 ended before R2"
                );
            }

            /*
             * R2 ended too early.
             */
            (Some(_), None) => {
                bail!(
                    "R2 ended before R1"
                );
            }

            /*
             * Both reads available.
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
                 * Make sure R1 and R2 still belong to
                 * the same pair.
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
                 * IMPORTANT:
                 *
                 * Both reads are classified against the
                 * SAME biological flank definitions.
                 *
                 * We search both orientations.
                 *
                 * Therefore we do NOT assume:
                 *
                 *     R1 = START -> END
                 *     R2 = END -> START
                 *
                 * and we do NOT assume that the input files
                 * are necessarily in the expected orientation.
                 */
                let r1_result =
                    classify_read(
                        &r1_seq,
                        start_flanks,
                        end_flanks,
                        args.max_hamming,
                    );

                let r2_result =
                    classify_read(
                        &r2_seq,
                        start_flanks,
                        end_flanks,
                        args.max_hamming,
                    );

                /*
                 * Normalize pair category.
                 *
                 * R1/R2 order does not matter here.
                 */
                let pair_category =
                    normalized_pair_category(
                        r1_result.category,
                        r2_result.category,
                    );

                pair_categories.add(
                    &r1_result,
                    &r2_result,
                );

                /*
                 * Keep physical R1/R2 lengths separate.
                 *
                 * This means the histogram still tells us
                 * which length belonged to the physical R1
                 * and which belonged to physical R2.
                 */
                histogram.add(
                    &pair_category,
                    r1_result.observed_length,
                    r2_result.observed_length,
                );

                /*
                 * ------------------------------------------------
                 * TSV
                 * ------------------------------------------------
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
                    r1_result.start_orientation,
                )?;

                write_optional_match(
                    &mut writer,
                    r1_result.end.as_ref(),
                    r1_result.end_orientation,
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
                    r2_result.start_orientation,
                )?;

                write_optional_match(
                    &mut writer,
                    r2_result.end.as_ref(),
                    r2_result.end_orientation,
                )?;

                writer.write_field(
                    r2_result
                        .observed_length
                        .to_string(),
                )?;

                /*
                 * Finish exactly one TSV row.
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
     * Write normalized pair categories.
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


/// Write normalized joint pair categories.
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
        "Reading biological flank definitions..."
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
        "Loaded {} biological START flanks and {} biological END flanks",
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
     * There are intentionally NO separate R1/R2 flank
     * definitions anymore.
     *
     * Every read is compared against:
     *
     *     biological START
     *     biological END
     *
     * in both:
     *
     *     forward
     *     reverse-complement
     *
     * orientations.
     *
     * Consequently:
     *
     *     START -> END
     *
     * and:
     *
     *     END -> START
     *
     * are both recognized as a complete "both" read.
     */
    process(
        &args,
        &start_flanks,
        &end_flanks,
    )
}