use super::permutation::PoseidonPermutation;
use crate::zq::Zq;

/// Errors that can occur while manipulating the underlying Poseidon state during
/// **absorption**.
#[derive(Debug)]
pub enum PoseidonError {
    /// Attempted to write to an index outside the state array.
    OutOfBounds(usize),
    /// Called with an *empty* input vector when one was expected.
    EmptyInput,
}

/// Errors returned by [`PoseidonSponge::squeeze`] when the caller’s request
/// is incompatible with the sponge *rate*.
#[derive(Debug, thiserror::Error)]
pub enum SpongeError {
    /// The requested number of field elements exceeds the sponge rate.
    #[error("Squeeze request of {requested} exceeds rate {rate}")]
    OversizedSqueeze { requested: usize, rate: usize },
}

/// Generic Poseidon **sponge construction**.
///
/// The sponge is parameterised over
/// * `OUTPUT_LENGTH` – number of field elements returned by [`compute_hash`].
/// * `RATE` – how many of the `WIDTH` state words are exposed during
///   *absorb*/*squeeze* (i.e. *capacity* = `WIDTH - RATE`).
/// * `WIDTH`, `ROUNDS`, `PARTIAL_ROUNDS`, `ALPHA` – forwarded to
///   [`PoseidonPermutation`].
///
/// **Usage pattern**:
/// 1. Create with [`new`], providing the message as `Vec<Zq>` and a
///    permutation instance.
/// 2. Call [`compute_hash`] to obtain the digest.
pub struct PoseidonSponge<
    const OUTPUT_LENGTH: usize,
    const RATE: usize,
    const WIDTH: usize,
    const ROUNDS: usize,
    const PARTIAL_ROUNDS: usize,
    const ALPHA: u64,
> {
    /// Buffer holding the message to be absorbed.
    input: Vec<Zq>,
    /// Output of length `OUTPUT_LENGTH`; initialised to all‑zero.
    output: [Zq; OUTPUT_LENGTH],
    /// Internal Poseidon permutation instance.
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
    /// Constructs a new sponge with the provided *message* and *permutation
    /// parameters*.
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

    /// **Absorb** the message into the sponge state according to the Poseidon
    /// rate/capacity split, forwarding to [`PoseidonPermutation::permute`] each
    /// time the rate section is filled. Returns an error only in pathological
    /// circumstances (e.g. writing out of bounds if `RATE`/`WIDTH` are
    /// miss‑configured).
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

    /// **Squeeze** `OUTPUT_LENGTH` field elements from the sponge. This is
    /// cached so that subsequent calls return the same slice without
    /// recomputation.
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

    /// Computes (and caches) the Poseidon hash of the input message. Repeated
    /// calls are **free** after the first one.
    pub fn compute_hash(&mut self) -> &[Zq] {
        if self.output.iter().all(|&x| x == Zq::ZERO) {
            self.absorb().expect("absorb failed");
            self.output = self.squeeze().expect("squeeze failed");
        }
        &self.output
    }
}

#[cfg(test)]
mod tests {
    use super::{PoseidonPermutation, PoseidonSponge};
    use crate::zq::Zq;

    #[test]
    fn test_absorb_fills_rate_then_permute() {
        // Arrange
        // WIDTH = 4, RATE = 2 → capacity = 2
        let mut permutation = PoseidonPermutation::<4, 2, 0, 1>::new();
        permutation.mds = [[Zq::ZERO; 4]; 4];
        for i in 0..4 {
            permutation.mds[i][i] = Zq::ONE;
        }
        permutation.ark = [[Zq::ZERO; 4]; 2];
        let msg = vec![Zq::new(5), Zq::new(7)];
        let mut sponge = PoseidonSponge::<1, 2, 4, 2, 0, 1>::new(msg, permutation);

        // absorb is implicit inside compute_hash
        sponge.absorb().unwrap();

        // After absorption the last `RATE` lanes should hold the message as we
        // used an identity permutation.
        assert_eq!(sponge.permutation.state[2], Zq::new(5));
        assert_eq!(sponge.permutation.state[3], Zq::new(7));
    }
}
