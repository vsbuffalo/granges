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
    .map_joins(mean_score)?;
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
map_multiple  270.21 s         112.66 s                      58.3073
map_max       105.46 s         84.03 s                       20.3185
adjust        112.42 s         53.48 s                       52.4269
filter        114.23 s         77.96 s                       31.7512
map_min       116.22 s         78.93 s                       32.0839
flank         162.77 s         80.67 s                       50.4383
map_mean      107.05 s         81.89 s                       23.5073
map_sum       108.43 s         93.35 s                       13.9083
windows       408.03 s         72.58 s                       82.2121
map_median    108.57 s         87.32 s                       19.5731
 
# run 2
command       bedtools time    granges time      granges speedup (%)
------------  ---------------  --------------  ---------------------
map_multiple  293.24 s         103.66 s                      64.6495
map_max       117.84 s         82.39 s                       30.0855
adjust        110.09 s         51.63 s                       53.0999
filter        120.36 s         67.79 s                       43.6784
map_min       114.76 s         86.06 s                       25.0081
flank         160.20 s         75.69 s                       52.756
map_mean      116.97 s         85.12 s                       27.2331
map_sum       114.39 s         85.96 s                       24.8557
windows       418.87 s         65.13 s                       84.4515
map_median    112.35 s         82.81 s                       26.2995
```

Note however, there is a arithmetic scaling issue. Larger range datasets
(1,000,000) lead to map operations that are inefficient. This is almost surely 
due to `data/operations.rs`. Here is an example benchmark:

```
command       bedtools time    granges time      granges speedup (%)
------------  ---------------  --------------  ---------------------
map_multiple  26.89 min        688.60 s                     57.3189
map_max       21.34 min        21.54 min                    -0.94889  *
adjust        21.35 min        499.54 s                     60.9959
filter        27.01 min        20.38 min                    24.5583
map_min       935.69 s         17.48 min                   -12.1008   * 
flank         30.26 min        943.36 s                     48.0353
map_mean      20.13 min        931.21 s                     22.913
map_sum       17.99 min        18.44 min                    -2.5225   *
windows       507.84 s         84.79 s                      83.304
map_median    16.90 min        17.81 min                    -5.36811  *
```

So median, sum, min, and max are slower.
