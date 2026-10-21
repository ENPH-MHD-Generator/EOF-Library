//! include/mhd_plasma.h is generated from the C ABI, and must be current.

use std::path::Path;

#[test]
fn the_c_header_matches_the_library() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let config = cbindgen::Config::from_file(root.join("cbindgen.toml")).expect("cbindgen.toml");
    let mut generated = Vec::new();
    cbindgen::Builder::new()
        .with_crate(root)
        .with_config(config)
        .generate()
        .expect("cbindgen could not read the crate")
        .write(&mut generated);
    let generated = String::from_utf8(generated).unwrap();
    let path = root.join("include/mhd_plasma.h");
    if std::env::var_os("UPDATE_HEADER").is_some() {
        std::fs::write(&path, &generated).unwrap();
        return;
    }
    let committed = std::fs::read_to_string(&path).unwrap_or_default();
    assert!(
        committed == generated,
        "include/mhd_plasma.h is stale: run `UPDATE_HEADER=1 cargo test` and commit it"
    );
}
