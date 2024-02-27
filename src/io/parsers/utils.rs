use std::path::Path;

/// Get the *base* extension to help infer filetype, which ignores compression-related
/// extensions (`.gz` and `.bgz`).
pub fn get_base_extension<P: AsRef<Path>>(filepath: P) -> Option<String> {
    let path = filepath.as_ref();

    // get the filename and split by '.'
    let parts: Vec<&str> = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("")
        .split('.')
        .collect();

    let ignore_extensions = ["gz", "bgz"];

    let has_ignore_extension = parts
        .last()
        .map_or(false, |ext| ignore_extensions.contains(ext));

    if parts.len() > 2 && has_ignore_extension {
        // if it's .gz, we return the second to last token,
        // e.g. path/foo.bed.gz would return bed
        Some(parts[parts.len() - 2].to_string())
    } else if parts.len() > 1 {
        // there is no .gz - return the last token.
        Some(parts[parts.len() - 1].to_string())
    } else {
        // no extension found
        None
    }
}
