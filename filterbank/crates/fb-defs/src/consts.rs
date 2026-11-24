pub const K_DM_PULSAR: f64 = 1.0 / 2.41e-4;
pub const K_DM_PHYSICAL: f64 = 4148.806_423_9;
pub const SECONDS_PER_DAY: f64 = 86_400.0;
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

pub fn dispersion_delay(dm: f64, f_mhz: f64, f_ref_mhz: f64) -> f64 {
    K_DM_PULSAR * dm * (f_mhz.powi(-2) - f_ref_mhz.powi(-2))
}

pub fn intrachannel_smearing(dm: f64, f_mhz: f64, chan_bw_mhz: f64) -> f64 {
    2.0 * K_DM_PULSAR * dm * chan_bw_mhz.abs() / f_mhz.powi(3)
}
