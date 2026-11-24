mod args;
mod convert;

use std::path::Path;
use std::process::ExitCode;

const HELP: &str = "fb — pulsar and transient search toolkit

commands:
  header FILE            print the canonical header (JSON) and any data gaps
  convert ...            PSRFITS/filterbank -> SIGPROC filterbank
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
        "-h" | "--help" | "help" => {
            print!("{HELP}\nconvert usage:\n  {}\n", convert::USAGE);
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
