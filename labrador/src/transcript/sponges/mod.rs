use crate::ring::{rq::Rq, zq::Zq};

pub trait Sponge: Default {
    fn absorb_zq(&mut self, input: &[Zq]);
    fn absorb_rq(&mut self, input: &[Rq]);
    fn squeeze_zq(&mut self, output_length: usize) -> Vec<Zq>;
    fn squeeze_rq(&mut self, output_length: usize) -> Vec<Rq>;
    fn squeeze_bits(&mut self, bit_length: usize) -> Vec<bool>;
    fn squeeze_bytes(&mut self, byte_length: usize) -> Vec<u8>;
}

pub mod shake;
// pub mod poseidon;
