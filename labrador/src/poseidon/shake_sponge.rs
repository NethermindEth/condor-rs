
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::ring::{rq::Rq, zq::Zq};

pub trait Sponge {
    fn absorb(&mut self, input: &[Rq]);
    fn squeeze(&mut self, output_lenght: usize) -> Vec<Rq>;
}


pub struct ShakeSponge {
    hasher: Shake256,
}

impl Default for ShakeSponge {
    fn default() -> Self {
        Self {
            hasher: Shake256::default()
        }
    }
}

impl Sponge for ShakeSponge {
    fn absorb(&mut self, input: &[Rq]) {
        // Convert Rq ector to u8
        let mut u8_version_input: Vec<u8> = Vec::new();
        for rq in input {
            for coeff in rq.get_coefficients() {
                u8_version_input.extend_from_slice(&coeff.get_value().to_be_bytes());
            }
        }
        self.hasher.update(&u8_version_input);
    }

    fn squeeze(&mut self, output_lenght: usize) -> Vec<Rq> {
        let mut reader = self.hasher.clone().finalize_xof();
        let mut output_buffer = vec![u8::default(); output_lenght * Rq::DEGREE * 4];
        reader.read(&mut output_buffer);


        let zq_values: Vec<Zq> = output_buffer
            .chunks_exact(4)
            .map(|chunk| u32::from_le_bytes(chunk.try_into().expect("Could not convert 4 u8 to one u32")))
            .map(|u32_value| Zq::new(u32_value))
            .collect();

        let result: Vec<Rq> = zq_values
            .chunks_exact(Rq::DEGREE)
            .map(|chunk| {
                let rq_input: [Zq; Rq::DEGREE] = chunk.try_into().expect("Chunk size is Rq::DEGREE");
                Rq::new(rq_input)
            })
            .collect();

        self.absorb(&result);
        result
    }
}
