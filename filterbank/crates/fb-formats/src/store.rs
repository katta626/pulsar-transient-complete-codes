use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use fb_defs::{DmPlan, ObsHeader, Result, RfiMask};
use serde::{de::DeserializeOwned, Serialize};

pub fn save_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    serde_json::to_writer_pretty(BufWriter::new(File::create(path)?), value)?;
    Ok(())
}

pub fn load_json<T: DeserializeOwned>(path: &Path) -> Result<T> {
    Ok(serde_json::from_reader(BufReader::new(File::open(path)?))?)
}

pub fn save_mask(path: &Path, m: &RfiMask) -> Result<()> {
    save_json(path, m)
}

pub fn load_mask(path: &Path) -> Result<RfiMask> {
    load_json(path)
}

pub fn save_header(path: &Path, h: &ObsHeader) -> Result<()> {
    save_json(path, h)
}

pub fn save_plan(path: &Path, p: &DmPlan) -> Result<()> {
    save_json(path, p)
}
