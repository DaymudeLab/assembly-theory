# assembly-theory

`assembly-theory` is an open-source, high-performance library for computing *assembly indices* of molecular structures (see, e.g., [Sharma et al., 2023](https://doi.org/10.1038/s41586-023-06600-9); [Walker et al., 2024](https://doi.org/10.1098/rsif.2024.0367)).
This crate is specific to the Rust library; see the
[GitHub repository](https://github.com/DaymudeLab/assembly-theory) for ways
to use `assembly-theory` as a standalone executable or as a Python package.


## Usage

Install the crate as usual:

```shell
cargo add assembly-theory
```

Load a molecule from a `.mol` file and calculate its assembly index:

```rust
use assembly_theory::{
    assembly::index,
    loader::parse_molfile_str
};

// Load a molecule from a `.mol` file.
let path = PathBuf::from(format!("anthracene.mol"));
let molfile = fs::read_to_string(path)?;
let anthracene = parse_molfile_str(&molfile).expect("Parsing failure.");

// Compute the molecule's assembly index using an efficient default algorithm.
assert_eq!(index(&anthracene), 6);
```

See [the documentation](https://docs.rs/assembly-theory) for a complete list of functions and usage exmaples.


## Citation

If you use this crate in your own scientific work, please consider citing us:

```bibtex
Coming soon!
```


## License

`assembly-theory` is licensed under the [Apache License, Version 2.0](https://choosealicense.com/licenses/apache-2.0/) or the [MIT License](https://choosealicense.com/licenses/mit/), at your option.

Unless you explicitly state otherwise, any contribution you intentionally submit for inclusion in this repository (as defined by Apache-2.0) shall be dual-licensed as above, without any additional terms or conditions.
