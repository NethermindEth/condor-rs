use crate::poseidon::PoseidonSponge;
use crate::poseidon::{PoseidonError, SpongeError};
use crate::zq::Zq;
use blake2::Blake2b256;
use blake2::Digest;
use rand::{random, rng, Rng};
use std::convert::TryInto;

pub trait Transcript {
    fn new() -> Self;
    fn absorb(&mut self, value: Zq) -> Result<(), PoseidonError>;
    fn get_challenge(&mut self) -> Result<Zq, SpongeError>;
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

    for (i, &x_val) in x_vals.iter().enumerate() {
        for (j, &y_val) in y_vals.iter().enumerate() {
            let denom = x_val + y_val;
            matrix[i][j] = denom.inv();
        }
    }

    matrix
}
// Generates ARK matrix
fn generate_secure_round_constants(rounds: usize, state_size: usize) -> Vec<Vec<Zq>> {
    // Create a matrix to store constants for each round and state position
    let mut ark = vec![vec![Zq::new(0); state_size]; rounds];

    // Iterate over the rows (rounds)
    for row in ark.iter_mut().take(rounds) {
        // Iterate over the states in the current row using `iter_mut` for mutable access
        for state in row.iter_mut().take(state_size) {
            // 4 random bytes to use as input (seed)
            let seed: [u8; 4] = random();

            // Hasher based on Blake2b
            let mut hasher = Blake2b256::new();
            hasher.update(seed);

            // Get the hash output
            let hash_result: [u8; 32] = hasher.finalize().into();

            // Take the first 4 bytes of the hash
            let hash_bytes: [u8; 4] = hash_result[..4]
                .try_into()
                .expect("can cast a slice into an array");

            // Convert those 4 bytes into a u32 value
            let constant = u32::from_le_bytes(hash_bytes);

            // Store the constant in the matrix
            *state = Zq::new(constant);
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
            generate_secure_round_constants(13, 16), // ark: Vec<Vec<Zq>>
        );

        Self { sponge }
    }

    fn absorb(&mut self, value: Zq) -> Result<(), PoseidonError> {
        self.sponge.absorb(&[value])?; // ? Propagates the error if it occurs
        Ok(())
    }

    fn get_challenge(&mut self) -> Result<Zq, SpongeError> {
        let squeezed_elements = self.sponge.squeeze(1);

        match squeezed_elements {
            Ok(vec) => Ok(vec[0] + Zq::ONE), // // Example: return first challenge + 1
            Err(err) => Err(err),            // Propagate the SpongeError
        }
    }
}
