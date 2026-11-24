pub mod bits;
pub mod fits;
pub mod psrfits;
pub mod reader;
pub mod sigproc;

pub use psrfits::{PsrfitsOptions, PsrfitsReader, ScalePolicy, ZeroOffPolicy};
pub use reader::{detect, open, open_with, Format, Gap, Reader};
pub use sigproc::{SigprocReader, SigprocWriter, WriteOptions};
