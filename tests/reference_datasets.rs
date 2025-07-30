//! Test assembly-theory correctness against all reference datasets.

use std::{collections::HashMap, ffi::OsStr, fs, path::Path};

use csv::Reader;

use assembly_theory::{
    assembly::{index_search, ParallelMode},
    bounds::Bound,
    canonize::CanonizeMode,
    enumerate::EnumerateMode,
    kernels::KernelMode,
    loader::parse_molfile_str,
    memoize::MemoizeMode,
};

fn load_ma_index(dataset: &str) -> HashMap<String, u32> {
    // Set up CSV reader for data/<dataset>/ma-index.csv.
    let ma_index_path = Path::new("data").join(dataset).join("ma-index.csv");
    let mut reader =
        Reader::from_path(ma_index_path).expect(&format!("{dataset}/ma-index.csv does not exist."));

    // Load assembly index records.
    let mut ma_index: HashMap<String, u32> = HashMap::new();
    for result in reader.records() {
        let record = result.expect("ma-index.csv is malformed.");
        let record = record.iter().collect::<Vec<_>>();
        ma_index.insert(
            record[0].to_string(),
            record[1]
                .to_string()
                .parse::<u32>()
                .expect("non-integer index"),
        );
    }

    // Return records.
    ma_index
}

fn test_reference_dataset(
    dataset: &str,
    enumerate_mode: EnumerateMode,
    canonize_mode: CanonizeMode,
    parallel_mode: ParallelMode,
    memoize_mode: MemoizeMode,
    kernel_mode: KernelMode,
    bounds: &[Bound],
) {
    // Load ground truth.
    let ma_index = load_ma_index(dataset);

    // Iterate over all .mol files in the dataset, computing the assembly index
    // of each one using the specified bounds. Track all molecules with indices
    // different than the ground truth.
    let mut incorrect_mols: Vec<(String, u32, u32)> = Vec::new();
    let mut paths: Vec<_> = fs::read_dir(Path::new("data").join(dataset))
        .unwrap()
        .filter_map(|r| r.ok())
        .collect();
    paths.sort_by_key(|p| p.path());
    for path in paths {
        // Only proceed if this is a .mol file.
        let name = path.path();
        if name.extension().and_then(OsStr::to_str) != Some("mol") {
            continue;
        }

        // Load the .mol file as an assembly_theory::molecule::Molecule.
        let mol = parse_molfile_str(
            &fs::read_to_string(name.clone()).expect(&format!("Could not read file {name:?}")),
        )
        .expect(&format!("Failed to parse {name:?}"));

        // Calculate the molecule's assembly index.
        let (index, _, _) = index_search(
            &mol,
            enumerate_mode,
            canonize_mode,
            parallel_mode,
            memoize_mode,
            kernel_mode,
            bounds,
        );

        // Compare calculated assembly index to ground truth.
        let molname = name.file_name().unwrap().to_str().unwrap().to_string();
        let true_index = ma_index[&molname];
        if index != true_index {
            incorrect_mols.push((molname, index, true_index));
        }
    }

    // If there are incorrect assembly indices, report and fail the test.
    let mut error_details = String::new();
    for (molname, index, true_index) in &incorrect_mols {
        error_details.push_str(&format!(
            "{molname}: assembly index {index} (assembly-theory) != {true_index} (ground truth)\n"
        ));
    }
    assert!(incorrect_mols.is_empty(), "{}", error_details);
}

/// Test enumeration modes individually without any search options.
///
/// This is very slow, so use gdb13_1201, the smallest dataset.
#[test]
fn enumeration() {
    for enumerate_mode in [EnumerateMode::Extend, EnumerateMode::GrowErode] {
        test_reference_dataset(
            "gdb13_1201",
            enumerate_mode,
            CanonizeMode::TreeNauty,
            ParallelMode::None,
            MemoizeMode::None,
            KernelMode::None,
            &[],
        );
    }
}

/// Test canonization modes individually without any search options.
///
/// This is very slow, so use gdb13_1201, the smallest dataset.
#[test]
fn canonization() {
    for canonize_mode in [CanonizeMode::Nauty, CanonizeMode::TreeNauty] {
        test_reference_dataset(
            "gdb13_1201",
            EnumerateMode::GrowErode,
            canonize_mode,
            ParallelMode::None,
            MemoizeMode::None,
            KernelMode::None,
            &[],
        );
    }
}

/// Test memoization modes individually.
///
/// Note that memoization with canonized keys can use any canonization mode.
/// Memoization does have some nontrivial interactions with parallelism, so we
/// test with parallelism. It does not use bounds, so use gdb13_1201 for speed.
///
/// TODO: It's possible that gdb13_1201 is a bad test for memoization since the
/// molecules may be so small that cached assembly states are never used again.
/// But coconut_55 molecules, which would be large enough for memoization to be
/// useful, are prohibitively slow without bounds.
#[test]
fn memoization() {
    for (memoize_mode, canonize_mode) in [
        (MemoizeMode::None, CanonizeMode::TreeNauty),
        (MemoizeMode::FragsIndex, CanonizeMode::TreeNauty),
        (MemoizeMode::CanonIndex, CanonizeMode::Nauty),
        (MemoizeMode::CanonIndex, CanonizeMode::TreeNauty),
    ] {
        test_reference_dataset(
            "gdb13_1201",
            EnumerateMode::GrowErode,
            canonize_mode,
            ParallelMode::DepthOne,
            memoize_mode,
            KernelMode::None,
            &[],
        );
    }
}

/// Test bounds individually.
///
/// The enumeration/canonization tests did not use bounds, so this tests one
/// actual bound at a time. Bounds don't have any race conditions, so this test
/// uses parallelism, but avoids any other search options.
#[test]
fn individual_bounds() {
    for bound in [
        Bound::Log,
        Bound::Int,
        Bound::VecSimple,
        Bound::VecSmallFrags,
    ] {
        test_reference_dataset(
            "checks",
            EnumerateMode::GrowErode,
            CanonizeMode::TreeNauty,
            ParallelMode::DepthOne,
            MemoizeMode::None,
            KernelMode::None,
            &[bound],
        );
    }
}

/// Test the application of all bounds simultaneously.
#[test]
fn all_bounds() {
    test_reference_dataset(
        "coconut_55",
        EnumerateMode::GrowErode,
        CanonizeMode::TreeNauty,
        ParallelMode::DepthOne,
        MemoizeMode::None,
        KernelMode::None,
        &[
            Bound::Log,
            Bound::Int,
            Bound::VecSimple,
            Bound::VecSmallFrags,
        ],
    );
}

/// Test the application of all bounds simultaneously with memoization.
#[test]
fn memoize_bounds() {
    test_reference_dataset(
        "coconut_55",
        EnumerateMode::GrowErode,
        CanonizeMode::TreeNauty,
        ParallelMode::DepthOne,
        MemoizeMode::CanonIndex,
        KernelMode::None,
        &[
            Bound::Log,
            Bound::Int,
            Bound::VecSimple,
            Bound::VecSmallFrags,
        ],
    );
}
