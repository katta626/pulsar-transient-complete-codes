#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BitOrder {
    MsbFirst,
    LsbFirst,
}

pub fn unpack_into(raw: &[u8], nbits: u32, order: BitOrder, signed: bool, out: &mut Vec<f32>) {
    match nbits {
        8 => {
            if signed {
                out.extend(raw.iter().map(|b| *b as i8 as f32));
            } else {
                out.extend(raw.iter().map(|b| *b as f32));
            }
        }
        16 => {
            out.extend(raw.chunks_exact(2).map(|c| {
                let v = [c[0], c[1]];
                if signed {
                    i16::from_ne_bytes(v) as f32
                } else {
                    u16::from_ne_bytes(v) as f32
                }
            }));
        }
        32 => out.extend(raw.chunks_exact(4).map(|c| f32::from_ne_bytes([c[0], c[1], c[2], c[3]]))),
        1 | 2 | 4 => {
            let per = 8 / nbits as usize;
            let mask = (1u8 << nbits) - 1;
            for &b in raw {
                for k in 0..per {
                    let shift = match order {
                        BitOrder::MsbFirst => 8 - nbits as usize * (k + 1),
                        BitOrder::LsbFirst => nbits as usize * k,
                    };
                    out.push(((b >> shift) & mask) as f32);
                }
            }
        }
        _ => {}
    }
}

pub fn unpack(raw: &[u8], nbits: u32, order: BitOrder, signed: bool) -> Vec<f32> {
    let mut out = Vec::with_capacity(raw.len() * 8 / nbits.max(1) as usize);
    unpack_into(raw, nbits, order, signed, &mut out);
    out
}

pub fn unpack_be16(raw: &[u8], signed: bool) -> Vec<f32> {
    raw.chunks_exact(2)
        .map(|c| if signed { i16::from_be_bytes([c[0], c[1]]) as f32 } else { u16::from_be_bytes([c[0], c[1]]) as f32 })
        .collect()
}

pub fn unpack_be32f(raw: &[u8]) -> Vec<f32> {
    raw.chunks_exact(4).map(|c| f32::from_be_bytes([c[0], c[1], c[2], c[3]])).collect()
}

pub fn pack(values: &[f32], nbits: u32, order: BitOrder, signed: bool) -> Vec<u8> {
    match nbits {
        8 => values
            .iter()
            .map(|v| if signed { v.round().clamp(-128.0, 127.0) as i8 as u8 } else { v.round().clamp(0.0, 255.0) as u8 })
            .collect(),
        16 => values
            .iter()
            .flat_map(|v| {
                if signed {
                    (v.round().clamp(-32768.0, 32767.0) as i16).to_le_bytes()
                } else {
                    (v.round().clamp(0.0, 65535.0) as u16).to_le_bytes()
                }
            })
            .collect(),
        32 => values.iter().flat_map(|v| v.to_le_bytes()).collect(),
        1 | 2 | 4 => {
            let per = 8 / nbits as usize;
            let maxv = ((1u32 << nbits) - 1) as f32;
            values
                .chunks(per)
                .map(|chunk| {
                    let mut b = 0u8;
                    for (k, v) in chunk.iter().enumerate() {
                        let q = v.round().clamp(0.0, maxv) as u8;
                        let shift = match order {
                            BitOrder::MsbFirst => 8 - nbits as usize * (k + 1),
                            BitOrder::LsbFirst => nbits as usize * k,
                        };
                        b |= q << shift;
                    }
                    b
                })
                .collect()
        }
        _ => Vec::new(),
    }
}
