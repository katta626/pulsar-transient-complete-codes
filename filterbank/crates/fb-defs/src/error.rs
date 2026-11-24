use std::fmt;

#[derive(Debug)]
pub enum FbError {
    Io(std::io::Error),
    Format(String),
    Incompatible(String),
    Unsupported(String),
    Invalid(String),
}

pub type Result<T> = std::result::Result<T, FbError>;

impl fmt::Display for FbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FbError::Io(e) => write!(f, "io error: {e}"),
            FbError::Format(s) => write!(f, "format error: {s}"),
            FbError::Incompatible(s) => write!(f, "incompatible inputs: {s}"),
            FbError::Unsupported(s) => write!(f, "unsupported: {s}"),
            FbError::Invalid(s) => write!(f, "invalid value: {s}"),
        }
    }
}

impl std::error::Error for FbError {}

impl From<std::io::Error> for FbError {
    fn from(e: std::io::Error) -> Self {
        FbError::Io(e)
    }
}

impl From<serde_json::Error> for FbError {
    fn from(e: serde_json::Error) -> Self {
        FbError::Format(e.to_string())
    }
}
