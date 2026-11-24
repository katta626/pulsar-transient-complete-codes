use serde::{Deserialize, Serialize};

use crate::header::ObsHeader;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Series {
    pub header: ObsHeader,
    pub dm: f64,
    pub f_ref: f64,
    pub data: Vec<f32>,
}

impl Series {
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn tsamp(&self) -> f64 {
        self.header.tsamp
    }

    pub fn mean_std(&self) -> (f64, f64) {
        let n = self.data.len().max(1) as f64;
        let m = self.data.iter().map(|v| *v as f64).sum::<f64>() / n;
        let v = self.data.iter().map(|x| (*x as f64 - m).powi(2)).sum::<f64>() / n;
        (m, v.sqrt())
    }
}
