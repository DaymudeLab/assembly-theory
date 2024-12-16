use crate::molecule;
use std::{
    io::{self, Error},
    path::PathBuf,
};

pub fn parse(p: &PathBuf) -> io::Result<molecule::Molecule> {
    Err(Error::new(
        io::ErrorKind::Unsupported,
        "Unimplemented operation!",
    ))
}
