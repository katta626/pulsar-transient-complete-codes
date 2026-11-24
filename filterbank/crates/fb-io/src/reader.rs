use std::path::Path;
use std::sync::Arc;

use fb_defs::{Block, FbError, ObsHeader, Result};

use crate::dada::DadaReader;
use crate::psrfits::{PsrfitsOptions, PsrfitsReader};
use crate::sigproc::SigprocReader;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Gap {
    pub start: u64,
    pub len: u64,
}

pub trait Reader: Send {
    fn header(&self) -> Arc<ObsHeader>;
    fn gaps(&self) -> &[Gap] {
        &[]
    }
    fn read(&mut self, start: u64, nsamp: usize) -> Result<Block>;

    fn blocks(&mut self, gulp: usize, overlap: usize) -> BlockIter<'_, Self>
    where
        Self: Sized,
    {
        BlockIter { reader: self, next: 0, gulp, overlap }
    }
}

pub struct BlockIter<'a, R: Reader> {
    reader: &'a mut R,
    next: u64,
    gulp: usize,
    overlap: usize,
}

impl<R: Reader> Iterator for BlockIter<'_, R> {
    type Item = Result<Block>;

    fn next(&mut self) -> Option<Self::Item> {
        let total = self.reader.header().nsamples;
        if self.next >= total {
            return None;
        }
        let b = self.reader.read(self.next, self.gulp + self.overlap);
        self.next += self.gulp as u64;
        Some(b)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Sigproc,
    Psrfits,
    Dada,
}

pub fn detect(path: &Path) -> Result<Format> {
    let mut f = std::fs::File::open(path)?;
    let mut buf = [0u8; 16];
    use std::io::Read;
    let n = f.read(&mut buf)?;
    if n >= 6 && &buf[..6] == b"SIMPLE" {
        return Ok(Format::Psrfits);
    }
    if n >= 16 && &buf[4..16] == b"HEADER_START" {
        return Ok(Format::Sigproc);
    }
    let text = String::from_utf8_lossy(&buf[..n]).to_ascii_uppercase();
    if text.starts_with("HDR_VERSION") || text.starts_with("HDR_SIZE") || path.extension().and_then(|e| e.to_str()) == Some("dada") {
        return Ok(Format::Dada);
    }
    Err(FbError::Unsupported(format!("cannot tell the format of {}", path.display())))
}

pub fn open<P: AsRef<Path>>(paths: &[P]) -> Result<Box<dyn Reader>> {
    open_with(paths, PsrfitsOptions::default())
}

pub fn open_with<P: AsRef<Path>>(paths: &[P], opts: PsrfitsOptions) -> Result<Box<dyn Reader>> {
    let first = paths.first().ok_or_else(|| FbError::Invalid("no input files".into()))?;
    match detect(first.as_ref())? {
        Format::Psrfits => Ok(Box::new(PsrfitsReader::open(paths, opts)?)),
        Format::Sigproc => {
            if paths.len() > 1 {
                return Err(FbError::Unsupported("multiple SIGPROC files: concatenate after reading each".into()));
            }
            Ok(Box::new(SigprocReader::open(first.as_ref())?))
        }
        Format::Dada => {
            if paths.len() > 1 {
                return Err(FbError::Unsupported("multiple DADA files: open them one at a time".into()));
            }
            Ok(Box::new(DadaReader::open(first.as_ref(), opts.sum_pols)?))
        }
    }
}
