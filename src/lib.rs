mod encode;
pub mod error;
mod index;
mod resonator;
mod weighted;

pub use bstr::{BStr, BString};
pub use error::ResonateError;
pub use resonator::HammingResonator;
pub use weighted::HammingResonatorWeighted;
