// Documentation

// Main Introduction

#![forbid(unsafe_code)]
#![deny(clippy::as_conversions)]
#![doc = include_str!("../../doc/mainpage-doc.md")]
// Arithmetic Circuit Translation
#![doc = include_str!("../../doc/arithmetic_circuit_translation.md")]
// Ajtai Commitment
#![doc = include_str!("../../doc/ajtai_commitment.md")]
// Projections
#![doc = include_str!("../../doc/projections.md")]
// Aggregation
#![doc = include_str!("../../doc/aggregation.md")]
// Amortization
#![doc = include_str!("../../doc/amortization.md")]

pub mod jl;
pub mod rq;
pub mod zq;

#[cfg(test)]
mod tests {

    #[test]
    fn test_for_workflow() {
        use base64::prelude::*;
        let _ = BASE64_STANDARD.encode([1, 2, 3, 4, 5]);
        assert_eq!(true, true)
    }
}
