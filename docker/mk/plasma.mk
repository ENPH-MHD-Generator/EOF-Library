# The mhd_plasma library (libs/solvers/MHDSolve/plasma, Rust): a static library
# with a C ABI that MHDSolve.so links. Rust's target glibc is 2.17, below this
# image's 2.23, and its objects are position-independent, so the archive links
# into the shared object as it is.
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

RUST_TOOLCHAIN := 1.92.0
export RUSTUP_HOME := /opt/rust/rustup
export CARGO_HOME  := /opt/rust/cargo
CARGO := $(CARGO_HOME)/bin/cargo

.PHONY: toolchain plasma

toolchain:
	wget -qO /tmp/rustup-init https://static.rust-lang.org/rustup/dist/x86_64-unknown-linux-gnu/rustup-init
	chmod +x /tmp/rustup-init
	# The base image sets HOME=/home/openfoam for a root build; rustup warns otherwise
	HOME=/root /tmp/rustup-init -y --no-modify-path --profile minimal --default-toolchain $(RUST_TOOLCHAIN)
	rm /tmp/rustup-init

plasma:
	cd $(EOF_SRC)/solvers/MHDSolve/plasma
	$(CARGO) build --release --locked
	install -D -m 0644 target/release/libmhd_plasma.a $(PLASMA_LIB)
