use crate::stats::quantile;

pub fn lags(radius: usize) -> Vec<isize> {
    let mut out = Vec::new();
    let mut k = 1usize;
    while k <= radius.max(1) {
        out.push(k as isize);
        out.push(-(k as isize));
        k = ((k as f64 * 1.5).floor() as usize).max(k + 1);
    }
    out
}

pub fn iqrm(x: &[f64], radius: usize, threshold: f64) -> Vec<bool> {
    let n = x.len();
    let mut cast = vec![0usize; n];
    let mut received = vec![0usize; n];
    let mut votes: Vec<(usize, usize)> = Vec::new();
    for k in lags(radius) {
        let d: Vec<f64> = (0..n)
            .map(|i| {
                let j = (i as isize - k).clamp(0, n as isize - 1) as usize;
                x[i] - x[j]
            })
            .collect();
        let mut tmp = d.clone();
        let q1 = quantile(&mut tmp, 0.25);
        let q2 = quantile(&mut tmp, 0.5);
        let q3 = quantile(&mut tmp, 0.75);
        let sigma = (q3 - q1) / 1.348_979_500_392_163_5;
        if sigma <= 0.0 {
            continue;
        }
        for i in 0..n {
            if d[i] - q2 > threshold * sigma {
                let j = (i as isize - k).clamp(0, n as isize - 1) as usize;
                if j != i {
                    votes.push((j, i));
                    cast[j] += 1;
                    received[i] += 1;
                }
            }
        }
    }
    let mut mask = vec![false; n];
    for (from, to) in votes {
        if cast[from] < received[to] {
            mask[to] = true;
        }
    }
    mask
}
