use criterion::{criterion_group, criterion_main, Criterion};
use csv::{self, ReaderBuilder};
use granges::io::parsers::Bed5Addition;
use granges::ranges::GenomicRangeRecord;
use granges::test_utilities::random_bed5file;
use granges::{prelude::*, Position};

#[derive(Debug, serde::Deserialize, PartialEq)]
struct Bed5Row {
    seqname: String,
    start: Position,
    end: Position,
    name: String,
    score: f64,
}

const BED_LENGTH: usize = 1_000_000;

fn bench_io_shootout(c: &mut Criterion) {
    // create the benchmark group
    let mut group = c.benchmark_group("adjust");

    // create the test data
    let input_bedfile = random_bed5file(BED_LENGTH);

    let genome = read_seqlens("tests_data/hg38_seqlens.tsv").unwrap();

    // configure the sample size for the group
    group.sample_size(10);

    // Bed5Iterator
    group.bench_function("bed5iterator", |b| {
        b.iter(|| {
            let iter = Bed5Iterator::new(input_bedfile.path()).unwrap();
            let gr = GRanges::from_iter(iter, &genome).unwrap();
            gr.len()
        });
    });

    // CSV
    group.bench_function("csv", |b| {
        b.iter(|| {
            let mut rdr = ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_path(input_bedfile.path())
                .unwrap();

            let iter = rdr.deserialize();
            let mut gr: GRanges<VecRangesIndexed, Vec<Bed5Addition>> = GRanges::new_vec(&genome);

            for result in iter {
                let row: GenomicRangeRecord<Bed5Addition> = result.unwrap();
                let data = Bed5Addition {
                    name: row.data.name,
                    score: row.data.score,
                };
                gr.push_range(&row.seqname, row.start, row.end, data)
                    .unwrap();
            }
        });
    });
}

criterion_group!(benches, bench_io_shootout,);
criterion_main!(benches);
