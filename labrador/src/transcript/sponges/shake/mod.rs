use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::ring::{rq::Rq, zq::Zq};
use crate::transcript::Sponge;

pub struct ShakeSponge {
    hasher: Shake256,
}

impl Sponge for ShakeSponge {
    fn default() -> Self {
        Self {
            hasher: Shake256::default(),
        }
    }
    fn absorb_zq(&mut self, input: &[Zq]) {
        // Convert Rq ector to u8
        let mut u8_version_input: Vec<u8> = Vec::new();
        for coeff in input {
            u8_version_input.extend_from_slice(&coeff.get_value().to_be_bytes());
        }
        self.hasher.update(&u8_version_input);
    }

    fn absorb_rq(&mut self, input: &[Rq]) {
        // Convert Rq ector to u8
        let mut u8_version_input: Vec<u8> = Vec::new();
        for rq in input {
            for coeff in rq.get_coefficients() {
                u8_version_input.extend_from_slice(&coeff.get_value().to_be_bytes());
            }
        }
        self.hasher.update(&u8_version_input);
    }

    fn squeeze_bits(&mut self, bit_length: usize) -> Vec<bool> {
        let byte_len = bit_length.div_ceil(8);
        let mut reader = self.hasher.clone().finalize_xof();
        let mut output_buffer = vec![u8::default(); byte_len];
        reader.read(&mut output_buffer);

        let mut result = Vec::with_capacity(bit_length);
        for byte in &output_buffer {
            let mut mask = 1u8;
            for _ in 0..8 {
                if result.len() == bit_length {
                    break;
                }
                result.push(byte & mask != 0);
                mask <<= 1;
            }
        }
        self.hasher.update(&output_buffer);
        result
    }

    fn squeeze_zq(&mut self, output_length: usize) -> Vec<Zq> {
        let mut reader = self.hasher.clone().finalize_xof();
        let mut output_buffer = vec![u8::default(); output_length * 4];
        reader.read(&mut output_buffer);

        let zq_values: Vec<Zq> = output_buffer
            .chunks_exact(4)
            .map(|chunk| {
                u32::from_le_bytes(chunk.try_into().expect("Could not convert 4 u8 to one u32"))
            })
            .map(Zq::new)
            .collect();

        self.absorb_zq(&zq_values);
        zq_values
    }

    fn squeeze_rq(&mut self, output_length: usize) -> Vec<Rq> {
        let mut reader = self.hasher.clone().finalize_xof();
        let mut output_buffer = vec![u8::default(); output_length * Rq::DEGREE * 4];
        reader.read(&mut output_buffer);

        let zq_values: Vec<Zq> = output_buffer
            .chunks_exact(4)
            .map(|chunk| {
                u32::from_le_bytes(chunk.try_into().expect("Could not convert 4 u8 to one u32"))
            })
            .map(Zq::new)
            .collect();

        let result: Vec<Rq> = zq_values
            .chunks_exact(Rq::DEGREE)
            .map(|chunk| {
                let rq_input: [Zq; Rq::DEGREE] =
                    chunk.try_into().expect("Chunk size is Rq::DEGREE");
                Rq::new(rq_input)
            })
            .collect();

        self.absorb_rq(&result);
        result
    }

    fn squeeze_bytes(&mut self, byte_length: usize) -> Vec<u8> {
        let mut reader = self.hasher.clone().finalize_xof();
        let mut output_buffer = vec![u8::default(); byte_length];
        reader.read(&mut output_buffer);
        self.hasher.update(&output_buffer);
        output_buffer
    }
}

#[cfg(test)]
mod test_sponge_correctness {
    use super::*;
    use crate::ring::zq::UniformZq;
    use rand::{distr::uniform::UniformSampler, rng};

    fn random_zq_vector(n: usize) -> Vec<Zq> {
        let uniform = UniformZq::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();
        (0..n).map(|_| uniform.sample(&mut rng())).collect()
    }

    #[test]
    fn test_zq_sponge_execution() {
        let mut sponge = ShakeSponge::default();
        let polynomial_1 = Rq::new([Zq::new(2); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(5); Rq::DEGREE]);
        let polynomial_3 = Rq::new([Zq::new(1); Rq::DEGREE]);

        sponge.absorb_rq(&[polynomial_1, polynomial_2, polynomial_3]);
        let result = sponge.squeeze_rq(1);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_rq_sponge_execution() {
        let mut sponge = ShakeSponge::default();
        let input: Vec<Zq> = random_zq_vector(64);
        sponge.absorb_zq(&input);
        let result = sponge.squeeze_rq(1);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_zq_squeeze_output_size() {
        let mut sponge = ShakeSponge::default();
        let input: Vec<Zq> = random_zq_vector(64);
        sponge.absorb_zq(&input);
        let result = sponge.squeeze_zq(1000);
        assert_eq!(result.len(), 1000);
    }

    #[test]
    fn test_rq_squeeze_output_size() {
        let mut sponge = ShakeSponge::default();
        let polynomial_1 = Rq::new([Zq::new(2); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(5); Rq::DEGREE]);
        let polynomial_3 = Rq::new([Zq::new(1); Rq::DEGREE]);

        sponge.absorb_rq(&[polynomial_1, polynomial_2, polynomial_3]);
        let result = sponge.squeeze_rq(1000);
        assert_eq!(result.len(), 1000);
    }

    #[test]
    fn test_absorb_zq_is_deterministic() {
        let input: Vec<Zq> = random_zq_vector(64);
        let mut s1 = ShakeSponge::default();
        let mut s2 = ShakeSponge::default();
        s1.absorb_zq(&input);
        s2.absorb_zq(&input);
        let out1 = s1.squeeze_zq(8);
        let out2 = s2.squeeze_zq(8);
        assert_eq!(out1, out2);
    }

    #[test]
    fn test_absorb_rq_is_deterministic() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);
        let polynomial_3 = Rq::new([Zq::new(9891741); Rq::DEGREE]);

        let mut sponge1 = ShakeSponge::default();
        sponge1.absorb_rq(&[
            polynomial_1.clone(),
            polynomial_2.clone(),
            polynomial_3.clone(),
        ]);

        let mut sponge2 = ShakeSponge::default();
        sponge2.absorb_rq(&[polynomial_1, polynomial_2, polynomial_3]);
        assert_eq!(sponge1.squeeze_rq(1), sponge2.squeeze_rq(1))
    }

    #[test]
    fn test_zq_successive_squeezes_is_unique() {
        let mut s = ShakeSponge::default();
        let input: Vec<Zq> = random_zq_vector(32);
        s.absorb_zq(&input);
        let o1 = s.squeeze_zq(8);
        let o2 = s.squeeze_zq(8);
        assert_ne!(o1, o2);
    }

    #[test]
    fn test_rq_successive_squeezes_is_unique() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);

        let mut sponge = ShakeSponge::default();
        sponge.absorb_rq(&[polynomial_1, polynomial_2]);
        let result1 = sponge.squeeze_rq(1);
        let result2 = sponge.squeeze_rq(1);
        assert_ne!(result1, result2)
    }

    #[test]
    fn test_zq_output_can_be_large() {
        let mut s = ShakeSponge::default();
        let input: Vec<Zq> = random_zq_vector(16);
        s.absorb_zq(&input);
        let out = s.squeeze_zq(1000);
        assert_eq!(out.len(), 1000);
    }

    #[test]
    fn test_rq_output_can_be_large() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);

        let mut sponge = ShakeSponge::default();
        sponge.absorb_rq(&[polynomial_1, polynomial_2]);

        assert_eq!(sponge.squeeze_rq(1000).len(), 1000);
    }

    #[test]
    fn test_squeeze_message_order() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);

        let mut sponge1 = ShakeSponge::default();
        sponge1.absorb_rq(&[polynomial_1.clone(), polynomial_2.clone()]);

        let mut sponge2 = ShakeSponge::default();
        sponge2.absorb_rq(&[polynomial_2, polynomial_1]);
        assert_ne!(sponge1.squeeze_rq(1), sponge2.squeeze_rq(1))
    }

    #[test]
    fn test_different_inputs_diff_outputs() {
        let input1: Vec<Zq> = random_zq_vector(64);
        let input2: Vec<Zq> = random_zq_vector(64);
        let mut s1 = ShakeSponge::default();
        let mut s2 = ShakeSponge::default();
        s1.absorb_zq(&input1);
        s2.absorb_zq(&input2);
        let o1 = s1.squeeze_zq(8);
        let o2 = s2.squeeze_zq(8);
        assert_ne!(o1, o2);
    }

    /// Edge case: zero‑length squeeze should not panic and must return empty vecs.
    #[test]
    fn zero_length_squeeze() {
        let mut s = ShakeSponge::default();
        assert!(s.squeeze_zq(0).is_empty());
        assert!(s.squeeze_rq(0).is_empty());
    }

    #[test]
    fn squeeze_bits_deterministic_and_length() {
        let mut s1 = ShakeSponge::default();
        let mut s2 = ShakeSponge::default();
        s1.absorb_zq(&[Zq::new(41)]);
        s2.absorb_zq(&[Zq::new(41)]);

        let b1 = s1.squeeze_bits(123);
        let b2 = s2.squeeze_bits(123);
        assert_eq!(b1, b2); // deterministic
        assert_eq!(b1.len(), 123); // exact length
    }
}

#[cfg(test)]
mod test_sponge_randomness {
    // We should test squeeze outputs are random-looking
}
