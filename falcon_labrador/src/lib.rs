// Documentation

// Main Introduction

#![forbid(unsafe_code)]
#![deny(clippy::as_conversions)]
#![doc = include_str!("../../doc/falcon_labrador_docs/mainpage-doc.md")]

///Example function
pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
