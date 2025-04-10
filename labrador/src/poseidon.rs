// Poseidon Sponge Implementation
//
// This file implements the Poseidon sponge construction, a cryptographic hash function (Based on Arkwork's implementation).
//
// Key components:
// - Zq: Finite field elements
// - PoseidonSponge: Struct holding sponge state and parameters.
//
// Sponge Operation:
// - 'rate' for input/output, 'capacity' for security.
// - 'permute' mixes the state using rounds.
//
// Functions:
// - new(): Initializes the sponge.
// - apply_ark(), apply_s_box(), apply_mds(): Round operations.
// - permute(): Applies Poseidon permutation.
// - absorb(): Absorbs input.
// - squeeze(): Squeezes output.
//

use crate::zq::Zq;
#[derive(Debug, thiserror::Error)]
pub enum SpongeError {
    #[error("Squeeze request of {requested} exceeds rate {rate}")]
    OversizedSqueeze { requested: usize, rate: usize },
}
#[derive(Debug)]
pub enum PoseidonError {
    OutOfBounds(usize),
    EmptyInput,
}


#[derive(Clone, Debug)]
/// Struct for the Poseidon sponge
pub struct PoseidonSponge {
    // Sponge parameters
    state: Vec<Zq>,  // Current state of the sponge
    rate: usize,     // Rate (number of elements to absorb/squeeze)
    capacity: usize, // Capacity (for cryptographic security)
    // Poseidon parameters
    full_rounds: usize, // Number of full rounds (S-box is applied to every element of the state)
    partial_rounds: usize, // Number of partial rounds (S-box is applied to first element of the state)
    alpha: u64,            // Exponent for the S-box
    mds: Vec<Vec<Zq>>,     // MDS matrix [[Zq; RATE + CAPACITY]; RATE + CAPACITY]
    ark: Vec<Vec<Zq>>,     // Round constants [[Zq; RATE + CAPACITY]; FULL_ROUNDS + PARTIAL_ROUNDS]
}

impl PoseidonSponge {
    // Initialize Poseidon sponge parameters
    pub fn new(
        rate: usize,
        capacity: usize,
        full_rounds: usize,
        partial_rounds: usize,
        alpha: u64,
        mds: Vec<Vec<Zq>>,
        ark: Vec<Vec<Zq>>,
    ) -> Self {
        // State starts with a vector of zeros
        let state = vec![Zq::ZERO; rate + capacity];
        Self {
            state,
            rate,
            capacity,
            full_rounds,
            partial_rounds,
            alpha,
            mds,
            ark,
        }
    }

    // Apply the round constants (ARK) (against algebraic attacks)
    fn apply_ark(&self, state: &mut [Zq], round_num: usize) {
        for (i, elem) in state.iter_mut().enumerate() {
            // adds diffusion (non-linearity)
            *elem += self.ark[round_num][i];
        }
    }

    // Apply the S-box (x -> x^alpha)
    fn apply_s_box(&self, state: &mut [Zq], is_full_round: bool) {
        if is_full_round {
            for elem in state {
                // Apply transformation to each element of state
                *elem = elem.pow(self.alpha);
            }
        } else {
            // Apply transformation only to the first element
            state[0] = state[0].pow(self.alpha);
        }
    }

    // Apply the MDS (Maximum Distance Separable) matrix
    fn apply_mds(&self, state: &mut [Zq]) {
        let mut new_state = state.to_vec();

        // Matrix multiplication with the state for diffusion
        for (i, _cur) in state.iter().enumerate() {
            let mut sum = Zq::ZERO;
            for (j, elem) in state.iter().enumerate() {
                // matrix multiplication
                sum += *elem * self.mds[i][j];
            }
            new_state[i] = sum;
        }

        state.copy_from_slice(&new_state);
    }

    // Apply a permutation
    fn permute(&mut self) {
        let full_rounds_half = self.full_rounds / 2;
        let mut state = self.state.clone();
        // full round ark/S-box/mds
        for i in 0..full_rounds_half {
            self.apply_ark(&mut state, i);
            self.apply_s_box(&mut state, true);
            self.apply_mds(&mut state);
        }
        // partial round ark/S-box/mds
        for i in full_rounds_half..(full_rounds_half + self.partial_rounds) {
            self.apply_ark(&mut state, i);
            self.apply_s_box(&mut state, false); // Apply the S-box to the first element only
            self.apply_mds(&mut state);
        }
        // full round ark/S-box/mds
        for i in (full_rounds_half + self.partial_rounds)..(self.full_rounds) {
            self.apply_ark(&mut state, i);
            self.apply_s_box(&mut state, true);
            self.apply_mds(&mut state);
        }

        self.state = state;
    }

    // Absorb input into the sponge
    pub fn absorb(&mut self, input: &[Zq]) -> Result<(), PoseidonError> {
        let mut remaining = input;
        let mut position = 0;
    
        while !remaining.is_empty() {
            let available_space = self.rate - position;
            let elements_to_absorb = remaining.len().min(available_space);
            let elements_to_process = &remaining[..elements_to_absorb];
            
            for (i, elem) in elements_to_process.iter().enumerate() {
                let target_index = self.capacity + position + i;
                if target_index >= self.state.len() {
                    return Err(PoseidonError::OutOfBounds(target_index));
                }
                self.state[target_index] += *elem;
            }
            
            remaining = &remaining[elements_to_absorb..];
            
            if remaining.is_empty() {
                self.permute();
                return Ok(());
            }
            
            self.permute();
            position = 0;
        }
        
        Ok(())
    }
    

    // Squeeze output from the sponge
    pub fn squeeze(&mut self, num_elements: usize) -> Result<Vec<Zq>, SpongeError> {
        // Check if the requested number of elements exceeds the rate.
        if num_elements > self.rate {
            return Err(SpongeError::OversizedSqueeze {
                requested: num_elements,
                rate: self.rate,
            });
        }

        let mut output = vec![Zq::ZERO; num_elements];
        let mut remaining_elements = output.as_mut_slice();
        let mut squeeze_index = 0;
        // Loop until we return the expected amount
        loop {
            if squeeze_index + remaining_elements.len() <= self.rate {
                remaining_elements.copy_from_slice(
                    &self.state[self.capacity + squeeze_index
                        ..self.capacity + remaining_elements.len() + squeeze_index],
                );
                return Ok(output);
            }
            // Reset squeeze index in case we need to squeeze again to fill the output
            let num_elements_squeezed = self.rate - squeeze_index;
            remaining_elements[..num_elements_squeezed].copy_from_slice(
                &self.state[self.capacity + squeeze_index
                    ..self.capacity + num_elements_squeezed + squeeze_index],
            );

            remaining_elements = &mut remaining_elements[num_elements_squeezed..];
            self.permute();
            squeeze_index = 0;
        }
    }
}
