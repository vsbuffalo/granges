use std::path::PathBuf;

use granges::io::io::read_seqlens;
use granges::prelude::*;



/// Adjust the genomic ranges in a bedfile by some specified amount.
// NOTE: we don't do build the full GRanges objects here, for efficiency.
// But it would be a good benchmark to see how much slower that is.
pub fn granges_adjust(bedfile: &PathBuf, seqlens: &PathBuf, both: &i32, output: Option<&PathBuf>) -> Result<(), GRangesError> {
    let genome = read_seqlens(seqlens)?;

    let bedlike_iterator = BedlikeIterator::new(bedfile)?;

    for record in bedlike_iterator {
        let range = record?;
        let seqname = range.seqname;
        let length = *genome.get(&seqname).ok_or(GRangesError::MissingSequence(seqname))?;
    }
    Ok(())
}
