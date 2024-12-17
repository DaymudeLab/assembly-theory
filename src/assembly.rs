use crate::molecule::Molecule;

fn top_down_search(m: &Molecule) {
    for (left, right) in m.partitions().unwrap() {
        println!("left: {:?}", left);
        println!("right: {:?}", right);
    }
}

// Compute the assembly index of a molecule
pub fn index(m: &Molecule) -> u32 {
    top_down_search(m);
    0
}
