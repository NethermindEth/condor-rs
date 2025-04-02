use crate::poseidon::PoseidonSponge;
use crate::zq::Zq;
use nalgebra::DMatrix;
use rand::{Rng, rngs::ThreadRng, thread_rng, random,rng}; 
use blake2::{Blake2b, Digest};
use generic_array::GenericArray; // Ensure this is imported for GenericArray
use typenum::U32;
use std::convert::TryInto;

pub trait Transcript {
    fn new() -> Self;
    fn absorb(&mut self, value: Zq);
    fn get_challenge(&mut self) -> Zq;
}

// Function to compute determinant mod q
fn determinant_mod_q(matrix: &Vec<Vec<Zq>>, q: u128) -> Zq {
    let size = matrix.len();
    let data: Vec<f64> = matrix.iter().flatten().map(|z| z.to_u128() as f64).collect();
    let mat = DMatrix::from_row_slice(size, size, &data);

    let det = mat.determinant().round() as i128 % q as i128;
    Zq::new(((det + q as i128) % q as i128) as u32)
}

// Function to check if a matrix is invertible in Zq
fn is_invertible(matrix: &Vec<Vec<Zq>>, q: u128) -> bool {
    determinant_mod_q(matrix, q) != Zq::new(0)
}

/// Generates an MDS (Maximum Distance Separable) matrix
fn mds_matrix(size: usize, q: u128) -> Vec<Vec<Zq>> {
    let mut rng = rng(); // Use `rng()` instead of `thread_rng()`
    let mut matrix: Vec<Vec<Zq>>;

    loop {
        // Generate a random matrix
        matrix = (0..size)
            .map(|_| (0..size).map(|_| Zq::new(rng.random_range(1..q) as u32)).collect())
            .collect();

        // Check if the matrix is invertible
        if is_invertible(&matrix, q) {
            break;
        }
    }

    matrix
}
fn generate_secure_round_constants(rounds: usize, state_size: usize) -> Vec<Vec<Zq>> {
    let mut ark = vec![vec![Zq::new(0); state_size]; rounds];

    // Set output size to 32 bytes (256 bits) for Blake2b
    for round in 0..rounds {
        for state in 0..state_size {
            let seed: [u8; 4] = random();  // Secure random 4-byte seed
            
            // Use Blake2b with the correct output size (32 bytes for 256-bit hash)
            let mut hasher = Blake2b::<U32>::new(); // Specify the output size type here
            hasher.update(seed);
            
            // Finalize the hash result
            let hash_result: GenericArray<u8, U32> = hasher.finalize();  // The hash result is the correct size
            
            // Extract the first 4 bytes of the hash result
            let hash_bytes: [u8; 4] = hash_result[..4].try_into().unwrap();  // Extract 4 bytes
            let constant = u32::from_le_bytes(hash_bytes);  // Convert to u32
            
            // Store the constant
            ark[round][state] = Zq::new(constant);
        }
    }

    ark
}
pub struct PoseidonTranscript {
    sponge: PoseidonSponge, // Transcrip will use the Poseidon sponge
}

impl Transcript for PoseidonTranscript {
    fn new() -> Self {
        // Define the Poseidon parameters (just examples for now)
        let sponge = crate::poseidon::PoseidonSponge::new(
            8,                                               // rate: usize
            8,                                               // capacity: usize
            5,                                               // full_rounds: usize
            8,                                               // partial_rounds: usize
            8,                                               // alpha: u64
            mds_matrix(16,2^(32) - 5),                      // mds: Vec<Vec<Zq>>
            generate_secure_round_constants(16, 16), // ark: Vec<Vec<Zq>>
        );

        //let sponge = PoseidonSponge::new(&params);  // Create the Poseidon sponge with params

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
