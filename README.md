![CI tests](https://github.com/vsbuffalo/granges/workflows/Rust/badge.svg)

## The GRanges Rust library and command line tool

GRanges is a Rust library for working with genomic ranges and their associated
data. It aims to make it easy to write extremely performant genomics tools that
work with genomic range data (e.g. BED, GTF/GFF, VCF, etc). Internally, GRanges
uses the *very* fast [coitrees](https://github.com/dcjones/coitrees/) interval
tree library written by Daniel C. Jones for overlap operations. In preliminary
benchmarks, GRanges tools can be 10%-30% faster than similar functionality in
[bedtools2](https://github.com/arq5x/bedtools2) (see benchmark and caveats
below). 

GRanges is inspired by ["tidy"](https://www.tidyverse.org) data analytics
workflows, as well as Bioconductor's
[GenomicRanges](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118)
and
[plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html).
GRanges uses a similar *method-chaining* pipeline approach to manipulate
genomic ranges, find overlapping genomic regions, and compute statistics.
For example, you could implement your own `bedtools map`-like functionality
in relatively few lines of code:

```rust
// Create the "right" GRanges object.
let right_gr = bed5_gr
    // Convert to interval trees.
    .into_coitrees()?
    // Extract out just the score from the additional BED5 columns.
    .map_data(|bed5_cols| {
        bed5_cols.score
    })?;

// Compute overlaps and combine scores into mean.
let results_gr = left_gr
    // Find overlaps
    .left_overlaps(&right_gr)?
    // Summarize overlap data
    .map_over_joins(mean_score)?;
```

where `mean_score()` is:

```rust
pub fn mean_score(join_data: CombinedJoinDataLeftEmpty<Option<f64>>) -> f64 {
    // Get the "right data" out of the join -- the BED5 scores.
    let overlap_scores: Vec<f64> = join_data.right_data.into_iter()
        // filter out missing values ('.' in BED)
        .filter_map(|x| x).collect();

    // Calculate the mean score.
    let score_sum: f64 = overlap_scores.iter().sum();
    score_sum / (overlap_scores.len() as f64)
}
```

Note that GRanges is a *compile-time* generic Rust library, so the code above
will be heavily optimized by the compiler. Rust uses *zero-cost abstractions*,
meaning high-level code like this is compiled and optimized so that it would be
just as performant as if it were written in a low-level language.

GRanges is *generic* in the sense that it works with *any* data container type
that stores data associated with genomic data: a `Vec<U>` of some type, an
[ndarray](https://docs.rs/ndarray/latest/ndarray/) `Array2`,
[polars](https://pola.rs) dataframe, etc. GRanges allows the user to write do
common genomics data processing tasks in a few lines of Rust, and then lets the
Rust compiler optimize it.

As a proof-of-concept, GRanges also provides the command line tool `granges`
built on this library's functionality. This command line tool is intended for
benchmarks against comparable command line tools and for large-scale
integration tests against other software to ensure that GRanges is bug-free.
The `granges` tool currently provides a subset of the features of other great
bioinformatics utilities like
[bedtools](https://bedtools.readthedocs.io/en/latest/). 

## Preliminary Benchmarks

In an attempt to combat "benchmark hype", this section details the results of
some preliminary benchmarks in an honest and transparent way. On our lab
server, here are two runs with 100,000 ranges per operation, and n = 100 samples:

```
# run 1
command       bedtools time    granges time      granges speedup (%)
------------  ---------------  --------------  ---------------------
map_multiple  293.41 s         129.99 s                     55.6963
map_max       131.97 s         111.39 s                     15.5938
adjust        127.36 s         59.25 s                      53.4811
filter        113.70 s         101.01 s                     11.1607
map_min       116.84 s         108.16 s                      7.42583
flank         150.95 s         86.25 s                      42.8602
map_mean      116.50 s         143.70 s                    -23.3524
map_sum       110.60 s         111.39 s                     -0.71377
windows       476.12 s         67.82 s                      85.7555
map_median    154.51 s         104.37 s                     32.4551

# run 2
command       bedtools time    granges time      granges speedup (%)
------------  ---------------  --------------  ---------------------
map_multiple  348.14 s         124.29 s                     64.2979
map_max       126.27 s         107.05 s                     15.2267
adjust        126.51 s         59.23 s                      53.1802
filter        114.66 s         97.28 s                      15.1542
map_min       118.61 s         113.02 s                      4.72037
flank         165.31 s         80.97 s                      51.0211
map_mean      133.39 s         108.22 s                     18.867
map_sum       120.28 s         128.68 s                     -6.98736
windows       476.22 s         90.89 s                      80.9151
map_median    106.83 s         104.37 s                      2.30808
```

The worse performance of `map_mean` in run 1, upon closer inspection, is
largely driven by [aberrant
replicates](https://github.com/vsbuffalo/granges/issues/2). (The server is
under heavy use by others, which adds variance to these benchmarks.) Note too
that unlike `bedtools`, `granges` *always* requires a genome file of chromosome
lengths for validation; this disk I/O adds to GRange's benchmark times. 

Here is another benchmark with smaller test datasets (10,000 ranges) and many
more replicates:

```
command       bedtools time    granges time      granges speedup (%)
------------  ---------------  --------------  ---------------------
map_multiple  133.41 s         77.81 s                      41.6736
adjust        60.23 s          40.06 s                      33.4794
map_min       67.11 s          59.31 s                      11.6281
map_mean      64.76 s          59.29 s                       8.43342
map_max       67.04 s          59.30 s                      11.5468
map_sum       66.89 s          60.01 s                      10.2835
map_median    65.63 s          59.47 s                       9.37911
flank         84.42 s          57.62 s                      31.7415
filter        77.95 s          56.65 s                      27.3318
windows       291.41 s         51.47 s                      82.3359
```

