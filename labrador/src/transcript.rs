use crate::poseidon::PoseidonSponge;
use crate::zq::Zq;
use blake2::{Blake2b, Digest};
use generic_array::GenericArray;
use rand::{random, rng, Rng};
use std::convert::TryInto;
use typenum::U32;

pub trait Transcript {
    fn new() -> Self;
    fn absorb(&mut self, value: Zq);
    fn get_challenge(&mut self) -> Zq;
}

/// Generates an MDS (Maximum Distance Separable) matrix using a Cauchy matrix
fn cauchy_mds_matrix(size: usize) -> Vec<Vec<Zq>> {
    let mut rng = rng();

    let mut x_vals = Vec::new();
    let mut y_vals = Vec::new();

    // Sample unique x values
    while x_vals.len() < size {
        let candidate = Zq::new(rng.random_range(1..u32::MAX));
        if !x_vals.contains(&candidate) {
            x_vals.push(candidate);
        }
    }

    // Sample unique y values with additional constraint
    while y_vals.len() < size {
        let candidate = Zq::new(rng.random_range(1..u32::MAX));
        let valid =
            !y_vals.contains(&candidate) && !x_vals.iter().any(|x| *x + candidate == Zq::ZERO);
        if valid {
            y_vals.push(candidate);
        }
    }

    // Construct the Cauchy matrix
    let mut matrix = vec![vec![Zq::ZERO; size]; size];
    for i in 0..size {
        for j in 0..size {
            let denom = x_vals[i] + y_vals[j];
            matrix[i][j] = denom.inv();
        }
    }

    matrix
}

fn generate_secure_round_constants(rounds: usize, state_size: usize) -> Vec<Vec<Zq>> {
    // Create a matrix to store constants for each round and state position
    let mut ark = vec![vec![Zq::new(0); state_size]; rounds];

    for round in 0..rounds {
        for state in 0..state_size {
            // 4 random bytes to use as input (seed)
            let seed: [u8; 4] = random();

            // Hasher based on Blake2b
            let mut hasher = Blake2b::<U32>::new();
            hasher.update(seed);

            // Get the hash output
            let hash_result: GenericArray<u8, U32> = hasher.finalize();

            // Take the first 4 bytes of the hash
            let hash_bytes: [u8; 4] = hash_result[..4].try_into().unwrap();

            // Convert those 4 bytes into a u32 value
            let constant = u32::from_le_bytes(hash_bytes);

            // Store the constant in the matrix as a Zq element
            ark[round][state] = Zq::new(constant);
        }
    }

    ark
}

pub struct PoseidonTranscript {
    sponge: PoseidonSponge,
}

impl Transcript for PoseidonTranscript {
    fn new() -> Self {
        // Parameters are just examples for now
        let sponge = crate::poseidon::PoseidonSponge::new(
            8,                                       // rate: usize
            8,                                       // capacity: usize
            5,                                       // full_rounds: usize
            8,                                       // partial_rounds: usize
            8,                                       // alpha: u64
            cauchy_mds_matrix(16),                   // mds: Vec<Vec<Zq>>
            generate_secure_round_constants(16, 16), // ark: Vec<Vec<Zq>>
        );

        Self { sponge }
    }

    fn absorb(&mut self, value: Zq) {
        self.sponge.absorb(&[value]); // Absorb the field value into the sponge
    }

    fn get_challenge(&mut self) -> Zq {
        // Squeeze the sponge to get the challenge (a Zq value)
        let squeezed_elements = self.sponge.squeeze(1); // we squeeze one number (as an example)
        squeezed_elements[0] + Zq::ONE // Return challenge + 1 as an example challenge
    }
}
