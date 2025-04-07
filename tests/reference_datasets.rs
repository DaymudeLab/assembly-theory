//! Test assembly-theory correctness against all reference datasets.

use csv::Reader;
use std::{
    collections::HashMap,
    ffi::OsStr,
    fs,
    path::Path
};

use assembly_theory::{
    loader,
    assembly::{
        index_search,
        serial_index_search,
        Bound
    }
};


fn load_ma_index(dataset: &str) -> HashMap<String, u32> {
    // Set up CSV reader for data/<dataset>/ma-index.csv.
    let ma_index_path = Path::new("data").join(dataset).join("ma-index.csv");
    let mut reader = Reader::from_path(ma_index_path)
        .expect(&format!("{dataset}/ma-index.csv does not exist."));
    
    // Load assembly index records.
    let mut ma_index: HashMap<String, u32> = HashMap::new();
    for result in reader.records() {
        let record = result.expect("ma-index.csv is malformed.");
        let record = record.iter().collect::<Vec<_>>();
        ma_index.insert(
            record[0].to_string(),
            record[1].to_string().parse::<u32>().expect("non-integer index"),
        );
    }

    // Return records.
    ma_index
}

fn test_reference_dataset(dataset: &str, bounds: &[Bound], serial: bool) {
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
        let mol = loader::parse_molfile_str(
            &fs::read_to_string(name.clone())
            .expect(&format!("Could not read file {name:?}"))
        ).expect(&format!("Failed to parse {name:?}"));

        // Calculate the molecule's assembly index.
        let (index, _, _) = if serial {
            serial_index_search(&mol, bounds)
        } else {
            index_search(&mol, bounds)
        };

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
        error_details.push_str(&format!("{molname}: assembly index {index} (assembly-theory) != {true_index} (ground truth)\n"));
    }
    assert!(incorrect_mols.is_empty(), "{}", error_details);
}

#[test]
fn gdb13_1201_naive() {
    test_reference_dataset("gdb13_1201", &vec![], false);
}

#[test]
fn gdb13_1201_logbound() {
    test_reference_dataset("gdb13_1201", &vec![Bound::Log], false);
}

#[test]
fn gdb13_1201_intbound() {
    test_reference_dataset("gdb13_1201", &vec![Bound::IntChain], false);
}

#[test]
fn gdb13_1201_allbounds() {
    let bounds = vec![Bound::IntChain,
                      Bound::VecChainSimple,
                      Bound::VecChainSmallFrags];
    test_reference_dataset("gdb13_1201", &bounds, false);
}

#[test]
fn gdb13_1201_naive_serial() {
    test_reference_dataset("gdb13_1201", &vec![], true);
}

#[test]
fn gdb13_1201_logbound_serial() {
    test_reference_dataset("gdb13_1201", &vec![Bound::Log], true);
}

#[test]
fn gdb13_1201_intbound_serial() {
    test_reference_dataset("gdb13_1201", &vec![Bound::IntChain], true);
}

#[test]
fn gdb13_1201_allbounds_serial() {
    let bounds = vec![Bound::IntChain,
                      Bound::VecChainSimple,
                      Bound::VecChainSmallFrags];
    test_reference_dataset("gdb13_1201", &bounds, true);
}

#[test]
#[ignore = "expensive test"]
fn gdb17_800_naive() {
    test_reference_dataset("gdb17_800", &vec![], false);
}

#[test]
#[ignore = "expensive test"]
fn gdb17_800_logbound() {
    test_reference_dataset("gdb17_800", &vec![Bound::Log], false);
}

#[test]
#[ignore = "expensive test"]
fn gdb17_800_intbound() {
    test_reference_dataset("gdb17_800", &vec![Bound::IntChain], false);
}

#[test]
#[ignore = "expensive test"]
fn gdb17_800_allbounds() {
    let bounds = vec![Bound::IntChain,
                      Bound::VecChainSimple,
                      Bound::VecChainSmallFrags];
    test_reference_dataset("gdb17_800", &bounds, false);
}

// #[test]
// #[ignore = "really expensive test"]
// fn checks_naive() {
//     test_reference_dataset("checks", &vec![], false);
// }

#[test]
fn checks_logbound() {
    test_reference_dataset("checks", &vec![Bound::Log], false);
}

#[test]
fn checks_intbound() {
    test_reference_dataset("checks", &vec![Bound::IntChain], false);
}

#[test]
fn checks_allbounds() {
    let bounds = vec![Bound::IntChain,
                      Bound::VecChainSimple,
                      Bound::VecChainSmallFrags];
    test_reference_dataset("checks", &bounds, false);
}

// #[test]
// #[ignore = "really expensive test"]
// fn coconut_220_naive() {
//     test_reference_dataset("coconut_220", &vec![], false);
// }

// #[test]
// #[ignore = "really expensive test"]
// fn coconut_220_logbound() {
//     test_reference_dataset("coconut_220", &vec![Bound::Log], false);
// }

// #[test]
// #[ignore = "really expensive test"]
// fn coconut_220_intbound() {
//     test_reference_dataset("coconut_220", &vec![Bound::IntChain], false);
// }

// #[test]
// #[ignore = "really expensive test"]
// fn coconut_220_allbounds() {
//     let bounds = vec![Bound::IntChain,
//                       Bound::VecChainSimple,
//                       Bound::VecChainSmallFrags];
//     test_reference_dataset("coconut_220", &bounds, false);
// }

