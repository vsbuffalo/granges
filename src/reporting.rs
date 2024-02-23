//! Types for standardized reports to the user about range operations.
//!
//! ⚠️: This is still under development.
//!
//! The goal of this is to encourage and facilitate developers of bioinformatics
//! tools to report information about potentially fragile operations, or inform
//! them of e.g. how many ranges were filtered out by some operation.
//!

/// The [`CommandOutput<U>`] type output is generic over some data output
/// from a command, and a [`Report`] that reports information to the user.
#[allow(unused)]
pub struct CommandOutput<U> {
    value: U,
    report: Report,
}

impl<U> CommandOutput<U> {
    pub fn new(value: U, report: Report) -> Self {
        Self { value, report }
    }
}

/// A type to (semi) standardize reporting to the user.
#[derive(Default)]
pub struct Report {
    entries: Vec<String>,
}

impl Report {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_issue(&mut self, message: String) {
        self.entries.push(message)
    }
}
