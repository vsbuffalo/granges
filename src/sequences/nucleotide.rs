//! Types and methods for working with per-basepair nucleotide sequence data.
//!
//! Currently this requires the [`noodles::fasta`] module, but their API is unstable
//! and may be a source of future pain.

use bytes::Bytes;
use genomap::GenomeMap;
use indexmap::IndexMap;
use noodles::core::Region;
use noodles::fasta::indexed_reader;
use noodles::fasta::{io::BufReadSeek, reader, record::Sequence, IndexedReader};
use std::ops::Deref;
use std::ops::Index;
use std::path::PathBuf;
use std::str;
use std::{cell::Ref, collections::HashSet, fmt};

use super::lazy::LazyLoader;
use crate::prelude::GRangesError;
use crate::ranges::try_range;
use crate::traits::Sequences;
use crate::{Position, INTERNAL_ERROR_MESSAGE};

/// A newtype around raw nucleotide [`Bytes`], for making it more
/// display and other operations more convenient.
#[derive(Clone, Debug, PartialEq)]
pub struct Nucleotides(Bytes);

impl fmt::Display for Nucleotides {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match str::from_utf8(&self.0) {
            Ok(s) => write!(f, "{}", s),
            Err(_) => Err(fmt::Error),
        }
    }
}

impl Deref for Nucleotides {
    type Target = Bytes;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<&Sequence> for Nucleotides {
    fn from(sequence: &Sequence) -> Self {
        let seq = Bytes::from(sequence.as_ref().to_vec());
        Nucleotides(seq)
    }
}

impl From<String> for Nucleotides {
    fn from(s: String) -> Self {
        let bytes = Bytes::from(s.into_bytes());
        Nucleotides(bytes)
    }
}

impl<'a> From<&'a str> for Nucleotides {
    fn from(s: &'a str) -> Self {
        let bytes = Bytes::from(s.as_bytes().to_vec());
        Nucleotides(bytes)
    }
}

/// A newtype around `Bytes` storage of nucleotide data
impl Nucleotides {
    /// Get the length of the nucleotide sequence.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Return whether this is an empty object.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return an [`NucleotideIterator`], which iterates through all
    /// bytes in a nucleotide sequence.
    pub fn iter(&self) -> NucleotideIterator<'_> {
        let inner = self.0.iter();
        NucleotideIterator { inner }
    }
}

/// Iterate over individual nucleotides.
pub struct NucleotideIterator<'a> {
    inner: std::slice::Iter<'a, u8>,
}

impl<'a> Iterator for NucleotideIterator<'a> {
    type Item = &'a u8;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}

/// [`NucleotideSequences`] for storing a whole genome's nucleotide sequence
/// data in-memory.
pub struct NucleotideSequences {
    data: GenomeMap<Nucleotides>,
}

impl NucleotideSequences {
    /// Load an entire FASTA file into memory, into a [`NucleotideSequences`] object.
    ///
    /// # Arguments
    /// * `filepath`: a path to the (possible gzipped) FASTA file.
    /// * `seqnames`: an optional subset of sequences to load.
    pub fn from_fasta(
        filepath: impl Into<PathBuf>,
        seqnames: Option<Vec<String>>,
    ) -> Result<Self, GRangesError> {
        let data = parse_fasta(filepath, seqnames)?;
        Ok(Self { data })
    }

    /// Retrieve an [`IndexMap`] of the sequence names and their lengths.
    pub fn seqlens(&self) -> IndexMap<String, Position> {
        self.data
            .iter()
            .map(|(k, v)| (k.clone(), v.len().try_into().unwrap()))
            .collect()
    }
}

/// Retrieve a nucleotide sequence by sequence name.
///
/// # Returns
/// Returns a reference to a `Nucleotides` object.
impl Index<&str> for NucleotideSequences {
    type Output = Nucleotides;

    fn index(&self, index: &str) -> &Self::Output {
        self.get_sequence(index).unwrap()
    }
}

impl Sequences for NucleotideSequences {
    type Container<'a> = &'a Nucleotides;
    type Slice<'a> = &'a [u8];

    /// Retrieve all sequence names.
    fn seqnames(&self) -> Vec<String> {
        self.seqlens().keys().cloned().collect()
    }

    /// Retrieve the [`Nucleotides`] for a particular sequence name.
    fn get_sequence(&self, seqname: &str) -> Result<Self::Container<'_>, GRangesError> {
        self.data
            .get(seqname)
            .ok_or(GRangesError::MissingSequenceName(seqname.to_string()))
    }

    /// Apply an arbitrary function to the specified region.
    ///
    /// # Arguments
    /// * `func`: a function that takes [`Bytes`] and processes them, returning a generic type `V`.
    /// * `seqname`: the sequence name of the region to apply the function to.
    /// * `start`: the start position of the region to apply the function to.
    /// * `end`: the end position of the region to apply the function to.
    fn region_map<V, F>(
        &self,
        func: &F,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<V, GRangesError>
    where
        F: Fn(Self::Slice<'_>) -> V,
    {
        let seq = self.get_sequence(seqname)?;
        let range = try_range(start, end, seq.len().try_into().unwrap())?;
        let data = &seq[range];
        Ok(func(data))
    }

    /// Get the length of a particular sequence.
    fn get_sequence_length(&self, seqname: &str) -> Result<Position, GRangesError> {
        let len: Position = self.get_sequence(seqname)?.len().try_into().unwrap();
        Ok(len)
    }
}

/// A lazy-loaded set of nucleotide sequences based on indexed FASTA files
#[derive(Debug)]
pub struct LazyNucleotideSequences {
    seqlens: IndexMap<String, Position>,
    lazy: LazyLoader<IndexedReader<Box<dyn BufReadSeek>>, Nucleotides, String>,
}

impl LazyNucleotideSequences {
    /// Create a new `LazyNucleotideSequences`, which can lazyily retrieve regions
    /// and whole sequences from an indexed FASTA file.
    ///
    /// # Arguments
    /// * `filepath` - the path to the bgzipped and indexed FASTA file.
    /// * `seqnames` - optional vector of sequences to consider (i.e. used for iterators over sequences).
    pub fn new(
        filepath: impl Into<PathBuf>,
        seqnames: Option<Vec<String>>,
    ) -> Result<Self, GRangesError> {
        let filepath = filepath.into();
        let reader = indexed_reader::Builder::default().build_from_path(filepath)?;
        let allowed_seqnames = seqnames.map(|c| c.into_iter().collect::<HashSet<_>>());

        let seqlens: IndexMap<String, Position> = reader
            .index()
            .iter()
            .filter_map(|r| {
                let name = String::from_utf8(r.name().to_vec()).unwrap_or_else(|_| {
                    panic!(
                        "{}\nError in converting FASTA name to a UTF8 String.",
                        &INTERNAL_ERROR_MESSAGE
                    )
                });
                if allowed_seqnames
                    .as_ref()
                    .map_or(true, |seqnames| seqnames.contains(&name))
                {
                    #[allow(clippy::useless_conversion)] // following clippy's suggestion causes
                    // compile time error
                    Some((name, r.length().try_into().unwrap()))
                } else {
                    None
                }
            })
            .collect();

        let lazy = LazyLoader::new(reader, move |reader, seqname: &String| {
            if allowed_seqnames
                .as_ref()
                .map_or(true, |seqnames| seqnames.contains(seqname))
            {
                // NOTE: noodles API changed to use &[u8] for unknown reasons,
                // which is rather annoying.
                let seqname_u8 = seqname.clone().into_bytes();
                let region = Region::new(seqname_u8, ..);
                let record = reader.query(&region)?;
                let seq = record.sequence();
                let nucs: Nucleotides = seq.into();
                Ok(nucs)
            } else {
                Err(GRangesError::MissingSequenceName(seqname.to_string()))
            }
        });

        Ok(Self { seqlens, lazy })
    }

    /// Check of the lazy-loading cache is empty.
    pub fn is_empty(&self) -> bool {
        self.lazy.is_empty()
    }

    /// Clear the lazy-loading cache.
    pub fn clear(&self) {
        self.lazy.clear()
    }

    /// Get an [`IndexMap`] of the sequence names and their lengths.
    pub fn seqlens(&self) -> IndexMap<String, Position> {
        self.seqlens.clone()
    }

    /// Return a `bool` indicating whether the specified `key` is cached.
    ///
    /// # Arguments
    /// * `key`: a `&str` key passed to the loader function.
    pub fn is_loaded(&self, key: &str) -> bool {
        self.lazy.is_loaded(&key.to_string())
    }
}

impl Sequences for LazyNucleotideSequences {
    type Container<'a> = Ref<'a, Nucleotides>;
    type Slice<'a> = &'a [u8];

    /// Retrieve all sequence names.
    fn seqnames(&self) -> Vec<String> {
        self.seqlens.keys().cloned().collect()
    }

    /// Retrieve the [`Nucleotides`] for a particular sequence name.
    fn get_sequence(&self, seqname: &str) -> Result<Self::Container<'_>, GRangesError> {
        self.lazy.get_data(&seqname.to_string())
    }

    /// Apply an arbitrary function to the specified region.
    ///
    /// # Arguments
    /// * `func`: a function that takes [`Bytes`] and processes them, returning a generic type `V`.
    /// * `seqname`: the sequence name of the region to apply the function to.
    /// * `start`: the start position of the region to apply the function to.
    /// * `end`: the end position of the region to apply the function to.
    fn region_map<V, F>(
        &self,
        func: &F,
        seqname: &str,
        start: Position,
        end: Position,
    ) -> Result<V, GRangesError>
    where
        F: for<'b> Fn(&'b [u8]) -> V,
    {
        let seq = self.get_sequence(seqname)?;
        let range = try_range(start, end, seq.len().try_into().unwrap())?;
        Ok(func(&seq[range]))
    }

    /// Get the length of a particular sequence.
    fn get_sequence_length(&self, seqname: &str) -> Result<Position, GRangesError> {
        self.seqlens
            .get(seqname)
            .ok_or(GRangesError::MissingSequenceName(seqname.to_string()))
            .copied()
    }
}

// Convert an `Option<Vec<String>>` into a `Option<HashSet<String>>`
fn option_vec_to_hashset(x: Option<Vec<String>>) -> Option<HashSet<String>> {
    x.map(HashSet::from_iter)
}

/// Use the [`noodles`] library to parse a FASTA file.
pub fn parse_fasta(
    filepath: impl Into<PathBuf>,
    seqnames: Option<Vec<String>>,
) -> Result<GenomeMap<Nucleotides>, GRangesError> {
    let seqnames_set = option_vec_to_hashset(seqnames);

    let filepath = filepath.into();

    let mut reader = reader::Builder.build_from_path(filepath)?;

    let mut sequences = GenomeMap::new();

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8(record.definition().name().to_vec())?;
        if seqnames_set
            .as_ref()
            .map_or(true, |keep_seqnames| keep_seqnames.contains(&name))
        {
            let seq = Bytes::from(record.sequence().as_ref().to_vec());
            sequences.insert(&name, Nucleotides(seq))?;
        }
    }

    Ok(sequences)
}

// some useful functions

/// Calculate the GC content of a byte slice, including IUPAC ambiguous codes.
///
/// # Arguments
/// * `seq` - a byte slice.
pub fn gc_content(seq: &[u8]) -> f64 {
    let gc_count = seq
        .iter()
        .filter(|&&base| {
            matches!(
                base.to_ascii_uppercase(),
                b'G' | b'C' | b'S' | b'V' | b'H' | b'D' | b'B' | b'N'
            )
        })
        .count() as f64;

    let total_count = seq.len() as f64;

    if total_count == 0.0 {
        0.0
    } else {
        gc_count / total_count
    }
}

/// Calculate the GC content of a byte slice, ignoring IUPAC ambiguous codes.
///
/// # Arguments
/// * `seq` - a byte slice.
pub fn gc_content_strict(seq: &[u8]) -> f64 {
    let gc_count = seq
        .iter()
        .filter(|&&base| matches!(base.to_ascii_uppercase(), b'G' | b'C'))
        .count() as f64;
    let total_count = seq
        .iter()
        .filter(|&&base| matches!(base.to_ascii_uppercase(), b'A' | b'T' | b'G' | b'C'))
        .count() as f64;

    if total_count == 0.0 {
        0.0
    } else {
        gc_count / total_count
    }
}

#[cfg(test)]
mod tests {
    use super::{gc_content_strict, LazyNucleotideSequences, NucleotideSequences};
    use crate::{granges::GRangesEmpty, sequences::nucleotide::Nucleotides, traits::Sequences};

    #[test]
    fn test_nucleotide_sequences() {
        let ref_file = "tests_data/sequences/test_case_01.fa.gz";
        let reference =
            NucleotideSequences::from_fasta(ref_file, None).expect("could not load reference");

        assert_eq!(
            *reference.get_sequence("chr1").unwrap(),
            Nucleotides::from("TTCACTACTATTAGTACTCACGGCGCAATA")
        );

        assert_eq!(reference.seqnames().len(), 2);
        assert_eq!(*reference.seqlens().get("chr1").unwrap(), 30);
    }

    #[test]
    fn test_lazyload() {
        let ref_file = "tests_data/sequences/test_case_01.fa.gz";
        let reference =
            LazyNucleotideSequences::new(ref_file, None).expect("could not load reference");

        assert_eq!(*reference.seqlens.get("chr1").unwrap(), 30);
        assert_eq!(*reference.seqlens.get("chr2").unwrap(), 100);

        let seqlens = reference.seqlens();

        // Test #1: nothing should be loaded yet
        assert_eq!(reference.is_loaded("chr1"), false);
        assert_eq!(reference.is_loaded("chr2"), false);

        // Test #2: validate the range summary by getting whole
        // chromosome sequence length through the apply funcs
        let seq1_len = seqlens.get("chr1").unwrap();

        fn get_len(seq: &[u8]) -> usize {
            seq.len()
        }

        let total_len = reference
            .region_map(&get_len, "chr1", 0, *seq1_len)
            .unwrap();
        assert_eq!(total_len, *seq1_len as usize);

        // Test #3: make sure the previous sequence is loaded still
        assert_eq!(reference.is_loaded("chr1"), true);

        // Test #4: make sure the next sequence is not loaded yet.
        assert_eq!(reference.is_loaded("chr2"), false);

        // Test #5: load the next sequence
        let chr2_len = seqlens.get("chr2").unwrap();
        let total_len = reference
            .region_map(&get_len, "chr2", 0, *chr2_len)
            .unwrap();
        assert_eq!(total_len, *chr2_len as usize);
        assert!(reference.is_loaded("chr2"));
    }

    #[test]
    fn lazy_region_apply_slice() {
        let ref_file = "tests_data/sequences/test_case_01.fa.gz";
        let reference =
            LazyNucleotideSequences::new(ref_file, None).expect("could not load reference");

        #[allow(non_snake_case)]
        // pycode for comparison
        // len([l for l in 'TTCACTACTATTAGTACTCACGGCGCAATA'[3:10] if l == 'C'])
        let total_Cs = reference
            .region_map(
                &|seq| seq.iter().filter(|c| **c == b'C').count(),
                "chr1",
                3,
                10,
            )
            .unwrap();
        assert_eq!(total_Cs, 2);

        #[allow(non_snake_case)]
        // pycode for comparison
        // len([l for l in 'TTCACTACTATTAGTACTCACGGCGCAATA'[3:10] if l == 'A'])
        let total_As = reference
            .region_map(
                &|seq| seq.iter().filter(|c| **c == b'A').count(),
                "chr1",
                3,
                10,
            )
            .unwrap();
        assert_eq!(total_As, 3);
    }

    #[test]
    fn test_map_into_granges() {
        let ref_file = "tests_data/sequences/test_case_01.fa.gz";
        let reference =
            LazyNucleotideSequences::new(ref_file, None).expect("could not load reference");

        let seqlens = reference.seqlens();

        let windows = GRangesEmpty::from_windows(&seqlens, 10, None, false).unwrap();

        let make_string = |seq: &[u8]| String::from_utf8(seq.to_vec()).unwrap();

        // our test: get subsequences (as Strings), with *no step*, and reconstruct the
        // original sequence
        let seqs_windows = reference
            .region_map_into_granges(&windows, &make_string)
            .unwrap();
        let subseqs = seqs_windows.data_refs_by_seqname().unwrap();
        let subseqs_chr1 = subseqs.get("chr1").unwrap();
        let chr1_reconstructed: Nucleotides = subseqs_chr1
            .iter()
            .map(AsRef::as_ref)
            .collect::<Vec<&str>>()
            .join("")
            .into();
        let chr1_seq = reference.get_sequence("chr1").unwrap();
        assert_eq!(chr1_reconstructed, *chr1_seq);
    }

    #[test]
    fn test_map_into_granges_gc() {
        let ref_file = "tests_data/sequences/test_case_01.fa.gz";
        let reference =
            LazyNucleotideSequences::new(ref_file, None).expect("could not load reference");

        let seqlens = reference.seqlens();

        let windows = GRangesEmpty::from_windows(&seqlens, 10, Some(5), false).unwrap();

        // our test: calc GC content
        let gc_windows = reference
            .region_map_into_granges(&windows, &gc_content_strict)
            .unwrap();
        let gc = gc_windows.data_by_seqname().unwrap();
        // TODO: double check? calc'd with help from py but should hand calc.
        assert_eq!(
            gc.get("chr1").unwrap().to_vec(),
            vec![0.3, 0.2, 0.3, 0.7, 0.6, 0.2]
        );
    }

    #[cfg(test)]
    #[cfg(feature = "human_tests")]
    fn test_lazyload_human() {
        use crate::{
            misc::assert_float_eq,
            prelude::*,
            sequences::nucleotide::{gc_content_strict, LazyNucleotideSequences},
        };

        let ref_file = "tests_data/human/hg38.analysisSet.fa.bgz";
        let file_exists = std::fs::metadata(ref_file).is_ok();

        let reference;
        if file_exists {
            let chroms: Vec<String> = (1..=22).map(|v| format!("chr{}", v)).collect();
            assert_eq!(chroms.len(), 22);
            reference = LazyNucleotideSequences::new(ref_file, Some(chroms))
                .expect("could not load reference");

            // expected hg38 lengths
            let expected_numseqs = reference.seqlens().unwrap().len();
            assert_eq!(expected_numseqs, 22);
        } else {
            println!(
                "human reference file '{}' not found, tests skipped.",
                ref_file
            );
            return;
        }

        fn get_len(seq: &[u8]) -> usize {
            seq.len()
        }

        let seqlens = reference.seqlens().unwrap();

        // Test #1: nothing should be loaded yet
        assert_eq!(reference.is_loaded("chr1"), false);

        // Test #2: validate the range summary by getting whole
        // chromosome sequence length through the apply funcs
        let chr1_len = seqlens.get("chr1").unwrap();
        // NOTE: the whole chromosome region is len(seq) - 1, since
        // it is 0-indexed, and *right-inclusive*
        let total_len = reference
            .region_apply(get_len, "chr1", 0, *chr1_len - 1)
            .unwrap();
        assert_eq!(total_len, *chr1_len as usize);

        // Test #3: make sure the previous sequence is loaded
        // still
        assert_eq!(reference.is_loaded("chr1"), true);

        // Test #4: make sure the next sequence is not loaded yet.
        assert_eq!(reference.is_loaded("chr2"), false);

        // Test #5: load the next sequence
        let chr2_len = seqlens.get("chr2").unwrap();
        let total_len = reference
            .region_apply(get_len, "chr2", 0, *chr2_len - 1)
            .unwrap();
        assert_eq!(total_len, *chr2_len as usize);

        // Test #6: hg specific -- test GC content
        // TODO should go in the human test module
        let gc = reference
            .region_apply(gc_content_strict, "chr1", 0, *chr1_len - 1)
            .unwrap();
        let expected_gc = 0.4172;
        assert_float_eq(gc, expected_gc, 0.0001);
    }
}
