pub fn median(v: &mut [f32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let n = v.len();
    let mid = n / 2;
    v.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let hi = v[mid];
    if n % 2 == 1 {
        hi
    } else {
        let lo = v[..mid].iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        0.5 * (lo + hi)
    }
}

pub fn quantile(v: &mut [f64], q: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let pos = q * (v.len() - 1) as f64;
    let i = pos.floor() as usize;
    let f = pos - i as f64;
    if i + 1 < v.len() {
        v[i] * (1.0 - f) + v[i + 1] * f
    } else {
        v[i]
    }
}

pub fn median_mad(v: &[f32]) -> (f32, f32) {
    let mut a = v.to_vec();
    let m = median(&mut a);
    for x in a.iter_mut() {
        *x = (*x - m).abs();
    }
    let mad = median(&mut a);
    (m, 1.4826 * mad)
}

#[derive(Debug, Clone, Default)]
pub struct ChannelStats {
    pub mean: Vec<f64>,
    pub std: Vec<f64>,
    pub skew: Vec<f64>,
    pub kurt: Vec<f64>,
    pub s1: Vec<f64>,
    pub s2: Vec<f64>,
}

pub fn channel_stats(data: &[f32], nsamp: usize, nchan: usize) -> ChannelStats {
    let mut s1 = vec![0f64; nchan];
    let mut s2 = vec![0f64; nchan];
    let mut s3 = vec![0f64; nchan];
    let mut s4 = vec![0f64; nchan];
    for t in 0..nsamp {
        let row = &data[t * nchan..(t + 1) * nchan];
        for c in 0..nchan {
            let x = row[c] as f64;
            let x2 = x * x;
            s1[c] += x;
            s2[c] += x2;
            s3[c] += x2 * x;
            s4[c] += x2 * x2;
        }
    }
    let n = nsamp.max(1) as f64;
    let mut st = ChannelStats { s1: s1.clone(), s2: s2.clone(), ..Default::default() };
    for c in 0..nchan {
        let m = s1[c] / n;
        let m2 = (s2[c] / n - m * m).max(0.0);
        let m3 = s3[c] / n - 3.0 * m * s2[c] / n + 2.0 * m.powi(3);
        let m4 = s4[c] / n - 4.0 * m * s3[c] / n + 6.0 * m * m * s2[c] / n - 3.0 * m.powi(4);
        let sd = m2.sqrt();
        st.mean.push(m);
        st.std.push(sd);
        st.skew.push(if sd > 0.0 { m3 / sd.powi(3) } else { 0.0 });
        st.kurt.push(if m2 > 0.0 { m4 / (m2 * m2) - 3.0 } else { 0.0 });
    }
    st
}
