//! Functionality for working with per-basepair data.
//!
//! The [`Sequences`] trait defines an abstraction over a per-basepair genomic range data, i.e.
//! when ranges do not need to be explicitly defined because they exhaustively cover the entire
//! genome (even if some values are "missing", e.g. a float NaN). The most common form of this type
//! of data is *nucleotide* sequence data (see [`nucleotide::NucleotideSequences`]) that defines
//! the reference genome nucleotide sequence. But, other type of sequence data exists, e.g. a
//! per-basepair numeric conservation score.
//!
//! ## Main Functionality
//!
//!  - Lazy-loading by sequence (e.g. chromosome) with the [`LazyLoader`] type.
//!  - If the `--features=ndarray` is set, one and two dimensional per-basepair
//!    numeric arrays/matrices *in-memory*, with [`NumericSequences1`] and
//!    [`NumericSequences2`], or *lazy loaded* with [`LazyNumericSequences2`]
//!    ([``LazyNumericSequences1`] will be added soon, please file a [GitHub
//!    issue](https://github.com/vsbuffalo/granges/issues) if you need it).
//!
//!
//! [`Sequences`]: crate::traits::Sequences
//! [`LazyLoader`]: crate::sequences::lazy
//!
//! [`NumericSequences1`]: crate::sequences::numeric::NumericSequences1
//! [`NumericSequences2`]: crate::sequences::numeric::NumericSequences2

pub mod lazy;
pub mod nucleotide;
#[cfg(feature = "ndarray")]
pub mod numeric;
