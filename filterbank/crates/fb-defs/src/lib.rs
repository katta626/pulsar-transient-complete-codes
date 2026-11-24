pub mod block;
pub mod candidates;
pub mod consts;
pub mod dmplan;
pub mod error;
pub mod header;
pub mod mask;
pub mod series;
pub mod time;

pub use block::Block;
pub use candidates::{PeriodicCandidate, SinglePulseCandidate};
pub use consts::*;
pub use dmplan::DmPlan;
pub use error::{FbError, Result};
pub use header::{DataOrder, FreqAxis, ObsHeader, PolType, Scaling};
pub use mask::RfiMask;
pub use series::Series;
pub use time::Epoch;

pub fn add_ext(base: &std::path::Path, ext: &str) -> std::path::PathBuf {
    let mut s = base.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    s.into()
}
