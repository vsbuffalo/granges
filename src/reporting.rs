//! Centralized reporting to the user about e.g. fragile operations.
//!

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
