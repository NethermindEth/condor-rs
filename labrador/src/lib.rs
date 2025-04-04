// Documentation

// Main Introduction
#![forbid(unsafe_code)]
#![deny(clippy::as_conversions)]
#![doc = include_str!("../doc/mainpage-doc.md")]
// Arithmetic Circuit Translation
#![doc = include_str!("../doc/arithmetic_circuit_translation.md")]
// Ajtai Commitment
#![doc = include_str!("../doc/ajtai_commitment.md")]
// Hierarchical Commitment
#![doc = include_str!("../doc/hierarchical_commitment.md")]
// Projections
#![doc = include_str!("../doc/projections.md")]
// Aggregation
#![doc = include_str!("../doc/aggregation.md")]
// Amortization
#![doc = include_str!("../doc/amortization.md")]

pub mod commitments;
pub mod core;
pub mod ring;
