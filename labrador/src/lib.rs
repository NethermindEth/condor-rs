// Documentation

// Main Introduction
#![forbid(unsafe_code)]
#![deny(clippy::as_conversions)]
#![doc = include_str!("../doc/labrador_docs/mainpage-doc.md")]
// Arithmetic Circuit Translation
#![doc = include_str!("../doc/labrador_docs/arithmetic_circuit_translation.md")]
// Ajtai Commitment
#![doc = include_str!("../doc/labrador_docs/ajtai_commitment.md")]
// Projections
#![doc = include_str!("../doc/labrador_docs/projections.md")]
// Aggregation
#![doc = include_str!("../doc/labrador_docs/aggregation.md")]
// Amortization
#![doc = include_str!("../doc/labrador_docs/amortization.md")]

pub mod rq;

pub mod zq;

pub mod ajtai_commitment;

pub mod rq_matrix;

pub mod rq_vector;

pub mod jl;
