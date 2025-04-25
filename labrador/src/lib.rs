// Documentation

// Main Introduction
#![forbid(unsafe_code)]
#![deny(clippy::as_conversions)]
#![doc = include_str!("../doc/mainpage-doc.md")]
// Arithmetic Circuit Translation
#![doc = include_str!("../doc/arithmetic_circuit_translation.md")]
// Ajtai Commitment
#![doc = include_str!("../doc/ajtai_commitment.md")]
// Projections
#![doc = include_str!("../doc/projections.md")]
// Aggregation
#![doc = include_str!("../doc/aggregation.md")]
// Amortization
#![doc = include_str!("../doc/amortization.md")]

pub mod rq;

pub mod zq;

pub mod ajtai_commitment;

pub mod rq_matrix;

pub mod rq_vector;

pub mod jl;

pub mod poseidon;

// pub mod transcript;
