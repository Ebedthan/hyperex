use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::fs;
use std::process::Command;

#[test]
fn file_do_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("hyperex")?;

    cmd.arg("tests/test.fa");
    cmd.assert().success();

    fs::remove_file("hyperex_out.fa")?;
    fs::remove_file("hyperex_out.gff")?;
    fs::remove_file("hyperex.log")?;

    Ok(())
}

#[test]
fn file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("hyperex")?;

    cmd.arg("test/file/doesnt/exists");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}
