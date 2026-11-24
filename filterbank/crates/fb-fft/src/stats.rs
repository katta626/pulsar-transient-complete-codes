pub fn ln_gamma(x: f64) -> f64 {
    let g = [
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = 0.999_999_999_999_809_9;
    let t = x + 7.5;
    for (i, gi) in g.iter().enumerate() {
        a += gi / (x + i as f64 + 1.0);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

pub fn ln_sum_power_sf(p: f64, m: u32) -> f64 {
    if p <= 0.0 {
        return 0.0;
    }
    let mut term = 0.0f64;
    let mut terms = Vec::with_capacity(m as usize);
    for j in 0..m {
        if j > 0 {
            term += p.ln() - (j as f64).ln();
        }
        terms.push(term);
    }
    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let s: f64 = terms.iter().map(|t| (t - mx).exp()).sum();
    (mx + s.ln() - p).min(0.0)
}

pub fn ln_erfc_sigma(ln_p: f64) -> f64 {
    if ln_p >= (0.5f64).ln() {
        return 0.0;
    }
    let mut x = (-2.0 * ln_p).sqrt().max(0.1);
    for _ in 0..80 {
        let lq = ln_gauss_sf(x);
        let pdf_ln = -0.5 * x * x - 0.5 * (2.0 * std::f64::consts::PI).ln();
        let dx = (lq - ln_p) * (lq - pdf_ln).exp();
        x += dx;
        if dx.abs() < 1e-10 {
            break;
        }
    }
    x
}

pub fn ln_gauss_sf(x: f64) -> f64 {
    if x < 5.0 {
        (0.5 * erfc(x / std::f64::consts::SQRT_2)).ln()
    } else {
        let x2 = x * x;
        -0.5 * x2 - x.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln() + (1.0 - 1.0 / x2 + 3.0 / (x2 * x2)).ln()
    }
}

pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t * (-z * z - 1.26551223
        + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
        .exp();
    if x >= 0.0 { r } else { 2.0 - r }
}

pub fn sigma_from_power(sum_power: f64, numharm: u32, ntrials: f64) -> f64 {
    let ln_p = ln_sum_power_sf(sum_power, numharm);
    let ln_pt = if ntrials > 1.0 {
        let p = ln_p.exp();
        if p * ntrials < 1e-6 { ln_p + ntrials.ln() } else { (1.0 - (1.0 - p).powf(ntrials)).max(1e-300).ln() }
    } else {
        ln_p
    };
    ln_erfc_sigma(ln_pt)
}
