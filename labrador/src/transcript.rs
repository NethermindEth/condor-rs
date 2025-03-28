// Example of transcript to try the Poseidon sponge
use crate::poseidon::PoseidonSponge;
use crate::zq::Zq;
pub trait Transcript {
    fn new() -> Self;
    fn absorb(&mut self, value: Zq);
    fn get_challenge(&mut self) -> Zq;
}

// MDS simple example (not cryptographically secure)
fn simple_mds_like_matrix(size: usize) -> Vec<Vec<Zq>> {
    let mut matrix = vec![vec![Zq::ZERO; size]; size];

    for i in 0..size {
        for j in 0..size {
            if i == j {
                matrix[i][j] = Zq::ONE;
            } else if i + j == size - 1 {
                matrix[i][j] = Zq::new(2);
            } else if i < j {
                matrix[i][j] = Zq::new(1);
            } else {
                matrix[i][j] = Zq::new(1);
            }
        }
    }
    matrix
}

// ARK simple example (not cryptographically secure)
fn generate_incremental_round_constants(
    rounds: usize,
    state_size: usize,
    start: u32,
) -> Vec<Vec<Zq>> {
    let mut ark = vec![vec![Zq::ZERO; state_size]; rounds];
    let mut counter = Zq::new(start);

    for i in 0..rounds {
        for j in 0..state_size {
            ark[i][j] = counter;
            counter += Zq::ONE;
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
            simple_mds_like_matrix(16),                      // mds: Vec<Vec<Zq>>
            generate_incremental_round_constants(16, 16, 1), // ark: Vec<Vec<Zq>>
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
