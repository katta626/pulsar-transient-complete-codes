use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use fb_defs::{FbError, Result};
use memmap2::Mmap;

const BLOCK: usize = 2880;
const CARD: usize = 80;

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Str(String),
    Int(i64),
    Float(f64),
    Bool(bool),
    None,
}

impl Value {
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Int(i) => Some(*i as f64),
            Value::Float(f) => Some(*f),
            Value::Str(s) => s.trim().parse().ok(),
            _ => None,
        }
    }

    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Int(i) => Some(*i),
            Value::Float(f) => Some(*f as i64),
            Value::Str(s) => s.trim().parse().ok(),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<String> {
        match self {
            Value::Str(s) => Some(s.trim().to_string()),
            Value::Int(i) => Some(i.to_string()),
            Value::Float(f) => Some(f.to_string()),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Header {
    pub cards: Vec<(String, Value)>,
    index: HashMap<String, usize>,
}

impl Header {
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.index.get(key).map(|&i| &self.cards[i].1)
    }

    pub fn f64(&self, key: &str) -> Option<f64> {
        self.get(key).and_then(|v| v.as_f64())
    }

    pub fn i64(&self, key: &str) -> Option<i64> {
        self.get(key).and_then(|v| v.as_i64())
    }

    pub fn str(&self, key: &str) -> Option<String> {
        self.get(key).and_then(|v| v.as_str()).filter(|s| !s.is_empty() && s != "*")
    }

    pub fn req_i64(&self, key: &str) -> Result<i64> {
        self.i64(key).ok_or_else(|| FbError::Format(format!("missing FITS keyword {key}")))
    }

    pub fn req_f64(&self, key: &str) -> Result<f64> {
        self.f64(key).ok_or_else(|| FbError::Format(format!("missing FITS keyword {key}")))
    }

    fn push(&mut self, key: String, value: Value) {
        if !self.index.contains_key(&key) {
            self.index.insert(key.clone(), self.cards.len());
        }
        self.cards.push((key, value));
    }
}

fn parse_value(raw: &str) -> Value {
    let s = raw.trim_start();
    if let Some(rest) = s.strip_prefix('\'') {
        let mut out = String::new();
        let mut chars = rest.chars().peekable();
        while let Some(c) = chars.next() {
            if c == '\'' {
                if chars.peek() == Some(&'\'') {
                    out.push('\'');
                    chars.next();
                } else {
                    break;
                }
            } else {
                out.push(c);
            }
        }
        return Value::Str(out.trim_end().to_string());
    }
    let v = s.split('/').next().unwrap_or("").trim();
    if v.is_empty() {
        return Value::None;
    }
    if v == "T" {
        return Value::Bool(true);
    }
    if v == "F" {
        return Value::Bool(false);
    }
    if let Ok(i) = v.parse::<i64>() {
        return Value::Int(i);
    }
    if let Ok(f) = v.replace(['D', 'd'], "E").parse::<f64>() {
        return Value::Float(f);
    }
    Value::Str(v.to_string())
}

#[derive(Debug, Clone)]
pub struct Column {
    pub name: String,
    pub repeat: usize,
    pub code: char,
    pub offset: usize,
    pub width: usize,
    pub dims: Option<Vec<usize>>,
}

impl Column {
    pub fn elem_size(code: char) -> usize {
        match code {
            'L' | 'B' | 'A' | 'X' => 1,
            'I' => 2,
            'J' | 'E' => 4,
            'K' | 'D' | 'C' | 'P' => 8,
            'M' | 'Q' => 16,
            _ => 1,
        }
    }
}

#[derive(Debug)]
pub struct Hdu {
    pub header: Header,
    pub data_offset: usize,
    pub data_len: usize,
}

impl Hdu {
    pub fn name(&self) -> String {
        self.header.str("EXTNAME").unwrap_or_else(|| "PRIMARY".into()).to_ascii_uppercase()
    }
}

pub struct FitsFile {
    pub mmap: Mmap,
    pub hdus: Vec<Hdu>,
}

impl FitsFile {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let mut hdus = Vec::new();
        let mut pos = 0usize;
        while pos + BLOCK <= mmap.len() {
            let (header, hlen) = read_header(&mmap[pos..])?;
            let data_offset = pos + hlen;
            let data_len = data_size(&header);
            hdus.push(Hdu { header, data_offset, data_len });
            pos = data_offset + data_len.div_ceil(BLOCK) * BLOCK;
        }
        if hdus.is_empty() {
            return Err(FbError::Format(format!("{} is not a FITS file", path.display())));
        }
        Ok(Self { mmap, hdus })
    }

    pub fn hdu(&self, name: &str) -> Option<&Hdu> {
        self.hdus.iter().find(|h| h.name() == name.to_ascii_uppercase())
    }

    pub fn primary(&self) -> &Header {
        &self.hdus[0].header
    }

    pub fn table(&self, name: &str) -> Result<Table<'_>> {
        let hdu = self.hdu(name).ok_or_else(|| FbError::Format(format!("no {name} HDU")))?;
        Table::new(hdu, &self.mmap)
    }
}

fn read_header(buf: &[u8]) -> Result<(Header, usize)> {
    let mut header = Header::default();
    let mut off = 0usize;
    loop {
        if off + CARD > buf.len() {
            return Err(FbError::Format("FITS header has no END card".into()));
        }
        let card = std::str::from_utf8(&buf[off..off + CARD]).unwrap_or("");
        off += CARD;
        let key = card[..8.min(card.len())].trim().to_string();
        if key == "END" {
            break;
        }
        if card.len() > 9 && &card[8..10] == "= " {
            header.push(key, parse_value(&card[10..]));
        } else if !key.is_empty() {
            header.push(key, Value::Str(card[8..].trim().to_string()));
        }
    }
    Ok((header, off.div_ceil(BLOCK) * BLOCK))
}

fn data_size(h: &Header) -> usize {
    let naxis = h.i64("NAXIS").unwrap_or(0) as usize;
    if naxis == 0 {
        return 0;
    }
    let bitpix = h.i64("BITPIX").unwrap_or(8).unsigned_abs() as usize;
    let mut n = 1usize;
    for i in 1..=naxis {
        n *= h.i64(&format!("NAXIS{i}")).unwrap_or(0) as usize;
    }
    let pcount = h.i64("PCOUNT").unwrap_or(0) as usize;
    let gcount = h.i64("GCOUNT").unwrap_or(1).max(1) as usize;
    bitpix / 8 * gcount * (pcount + n)
}

pub struct Table<'a> {
    pub header: &'a Header,
    pub row_len: usize,
    pub nrows: usize,
    pub columns: Vec<Column>,
    data: &'a [u8],
}

impl<'a> Table<'a> {
    fn new(hdu: &'a Hdu, mmap: &'a Mmap) -> Result<Self> {
        let h = &hdu.header;
        let row_len = h.req_i64("NAXIS1")? as usize;
        let nrows = h.req_i64("NAXIS2")? as usize;
        let ncol = h.req_i64("TFIELDS")? as usize;
        let mut columns = Vec::with_capacity(ncol);
        let mut offset = 0usize;
        for i in 1..=ncol {
            let name = h.str(&format!("TTYPE{i}")).unwrap_or_default().to_ascii_uppercase();
            let form = h.str(&format!("TFORM{i}")).ok_or_else(|| FbError::Format(format!("missing TFORM{i}")))?;
            let digits: String = form.chars().take_while(|c| c.is_ascii_digit()).collect();
            let repeat: usize = if digits.is_empty() { 1 } else { digits.parse().unwrap_or(1) };
            let code = form[digits.len()..].chars().next().unwrap_or('B');
            let width = if code == 'X' { repeat.div_ceil(8) } else { repeat * Column::elem_size(code) };
            let dims = h.str(&format!("TDIM{i}")).map(|d| {
                d.trim_matches(|c| c == '(' || c == ')').split(',').filter_map(|x| x.trim().parse().ok()).collect()
            });
            columns.push(Column { name, repeat, code, offset, width, dims });
            offset += width;
        }
        let end = hdu.data_offset + row_len * nrows;
        if end > mmap.len() {
            return Err(FbError::Format("binary table extends past end of file".into()));
        }
        Ok(Self { header: h, row_len, nrows, columns, data: &mmap[hdu.data_offset..end] })
    }

    pub fn column(&self, name: &str) -> Option<&Column> {
        self.columns.iter().find(|c| c.name == name)
    }

    pub fn cell(&self, row: usize, col: &Column) -> &'a [u8] {
        let start = row * self.row_len + col.offset;
        &self.data[start..start + col.width]
    }

    pub fn f64s(&self, row: usize, name: &str) -> Option<Vec<f64>> {
        let col = self.column(name)?;
        let b = self.cell(row, col);
        Some(match col.code {
            'D' => b.chunks_exact(8).map(|c| f64::from_be_bytes(c.try_into().unwrap())).collect(),
            'E' => b.chunks_exact(4).map(|c| f32::from_be_bytes(c.try_into().unwrap()) as f64).collect(),
            'J' => b.chunks_exact(4).map(|c| i32::from_be_bytes(c.try_into().unwrap()) as f64).collect(),
            'I' => b.chunks_exact(2).map(|c| i16::from_be_bytes(c.try_into().unwrap()) as f64).collect(),
            'K' => b.chunks_exact(8).map(|c| i64::from_be_bytes(c.try_into().unwrap()) as f64).collect(),
            'B' => b.iter().map(|v| *v as f64).collect(),
            _ => return None,
        })
    }

    pub fn f64_scalar(&self, row: usize, name: &str) -> Option<f64> {
        self.f64s(row, name).and_then(|v| v.first().copied())
    }
}
