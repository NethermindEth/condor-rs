use super::permutation::PoseidonPermutation;
use crate::zq::Zq;

#[derive(Debug)]
pub enum PoseidonError {
    OutOfBounds(usize),
    EmptyInput,
}

#[derive(Debug, thiserror::Error)]
pub enum SpongeError {
    #[error("Squeeze request of {requested} exceeds rate {rate}")]
    OversizedSqueeze { requested: usize, rate: usize },
}

pub struct PoseidonSponge<
    const OUTPUT_LENGTH: usize,
    const RATE: usize,
    const WIDTH: usize,
    const ROUNDS: usize,
    const PARTIAL_ROUNDS: usize,
    const ALPHA: u64,
> {
    input: Vec<Zq>,
    output: [Zq; OUTPUT_LENGTH],
    permutation: PoseidonPermutation<WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>,
}

impl<
        const OUTPUT_LENGTH: usize,
        const RATE: usize,
        const WIDTH: usize,
        const ROUNDS: usize,
        const PARTIAL_ROUNDS: usize,
        const ALPHA: u64,
    > PoseidonSponge<OUTPUT_LENGTH, RATE, WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>
{
    pub fn new(
        input: Vec<Zq>,
        permutation: PoseidonPermutation<WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>,
    ) -> Self {
        Self {
            input,
            output: [Zq::ZERO; OUTPUT_LENGTH],
            permutation,
        }
    }

    pub fn default(
        input: Vec<Zq>,
        permutation: PoseidonPermutation<WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>,
    ) -> Self {
        Self::new(input, permutation)
    }

    // Absorb input into the sponge
    fn absorb(&mut self) -> Result<(), PoseidonError> {
        let mut remaining = self.input.clone();
        let mut position = 0;

        while !remaining.is_empty() {
            let available_space = RATE - position;
            let elements_to_absorb = remaining.len().min(available_space);
            let elements_to_process = &remaining[..elements_to_absorb];

            for (i, elem) in elements_to_process.iter().enumerate() {
                let target_index = WIDTH - RATE + position + i;
                if target_index >= self.permutation.state.len() {
                    return Err(PoseidonError::OutOfBounds(target_index));
                }
                self.permutation.state[target_index] += *elem;
            }

            remaining = (remaining[elements_to_absorb..]).to_vec();

            if remaining.is_empty() {
                self.permutation.permute();
                return Ok(());
            }
            self.permutation.permute();
            position = 0;
        }

        Ok(())
    }

    // Squeeze output from the sponge
    fn squeeze(&mut self) -> Result<[Zq; OUTPUT_LENGTH], SpongeError> {
        // Check if the requested number of elements exceeds the rate.
        if OUTPUT_LENGTH > RATE {
            return Err(SpongeError::OversizedSqueeze {
                requested: OUTPUT_LENGTH,
                rate: RATE,
            });
        }

        let mut output = [Zq::ZERO; OUTPUT_LENGTH];
        let mut remaining_elements = output.as_mut_slice();
        let mut squeeze_index = 0;
        // Loop until we return the expected amount
        loop {
            if squeeze_index + remaining_elements.len() <= RATE {
                remaining_elements.copy_from_slice(
                    &self.permutation.state[WIDTH - RATE + squeeze_index
                        ..WIDTH - RATE + remaining_elements.len() + squeeze_index],
                );
                return Ok(output);
            }
            // Reset squeeze index in case we need to squeeze again to fill the output
            let num_elements_squeezed = RATE - squeeze_index;
            remaining_elements[..num_elements_squeezed].copy_from_slice(
                &self.permutation.state[WIDTH - RATE + squeeze_index
                    ..WIDTH - RATE + num_elements_squeezed + squeeze_index],
            );

            remaining_elements = &mut remaining_elements[num_elements_squeezed..];
            self.permutation.permute();
            squeeze_index = 0;
        }
    }

    pub fn compute_hash(&mut self) -> &[Zq] {
        if self.output.iter().all(|&x| x == Zq::ZERO) {
            self.absorb().expect("absorb failed");
            self.output = self.squeeze().expect("squeeze failed");
        }
        &self.output
    }

    pub fn reset(&mut self) {
        self.input = vec![];
        self.output = [Zq::ZERO; OUTPUT_LENGTH];
        self.permutation = PoseidonPermutation::new();
    }
}
