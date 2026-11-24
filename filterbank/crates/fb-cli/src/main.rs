mod args;
mod convert;
mod periodic;
mod search;

use std::path::Path;
use std::process::ExitCode;

const HELP: &str = "fb — pulsar and transient search toolkit

commands:
  header FILE            print the canonical header (JSON) and any data gaps
  convert ...            PSRFITS/filterbank -> SIGPROC filterbank
  rfi ...                build an RFI mask (JSON, optional PRESTO .mask)
  sp ...                 single-pulse search: RFI -> dedisperse -> boxcar -> cluster
  export-mask ...        convert masks between JSON and PRESTO .mask
  periodic ...           periodicity search: RFI -> dedisperse -> FFT/accel + FFA -> sift -> fold
  fft SERIES             FFT (+acceleration) search of one time series (.dat/.tim)
  ffa SERIES             fast folding algorithm search of one time series
  fold FILE              fold a series or filterbank at a given period
";

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = argv.first() else {
        eprint!("{HELP}");
        return ExitCode::from(2);
    };
    let rest = &argv[1..];
    let res = match cmd.as_str() {
        "header" => match rest.first() {
            Some(p) => convert::header(Path::new(p)),
            None => Err(fb_defs::FbError::Invalid("usage: fb header FILE".into())),
        },
        "convert" => convert::run(rest),
        "sp" => search::sp(rest),
        "rfi" => search::rfi(rest),
        "export-mask" => search::export_mask(rest),
        "periodic" => periodic::periodic(rest),
        "fft" => periodic::fft(rest),
        "ffa" => periodic::ffa(rest),
        "fold" => periodic::fold(rest),
        "-h" | "--help" | "help" => {
            print!(
                "{HELP}\nconvert:\n  {}\nsp:\n  {}\nrfi:\n  {}\nperiodic:\n  {}\nseries:\n  {}\n",
                convert::USAGE,
                search::SP_USAGE,
                search::RFI_USAGE,
                periodic::PERIODIC_USAGE,
                periodic::SERIES_USAGE
            );
            Ok(())
        }
        other => Err(fb_defs::FbError::Invalid(format!("unknown command '{other}'\n{HELP}"))),
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("fb: {e}");
            ExitCode::from(1)
        }
    }
}
