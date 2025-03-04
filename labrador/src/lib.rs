// Documentation

// Main Introduction

#![forbid(unsafe_code)]
#![deny(clippy::as_conversions)]
#![doc = include_str!("../../doc/labrador_docs/mainpage-doc.md")]
// Arithmetic Circuit Translation
#![doc = include_str!("../../doc/labrador_docs/arithmetic_circuit_translation.md")]
// Ajtai Commitment
#![doc = include_str!("../../doc/labrador_docs/ajtai_commitment.md")]
// Projections
#![doc = include_str!("../../doc/labrador_docs/projections.md")]
// Aggregation
#![doc = include_str!("../../doc/labrador_docs/aggregation.md")]
// Amortization
#![doc = include_str!("../../doc/labrador_docs/amortization.md")]

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
