use std::ops::Index;

use super::{
    constants::SpongeConstants,
    permutation::PoseidonPermutation,
};
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

pub struct PoseidonHash<S: SpongeConstants> {
    input: Vec<Zq>,
    output: Vec<Zq>,
    permutation: PoseidonPermutation<S>,
}

impl<S> PoseidonHash<S>
where
    S: SpongeConstants,
    S::ArkArray: Index<usize>,
    S::MdsArray: Index<usize>,
    <S::ArkArray as Index<usize>>::Output: Index<usize, Output = Zq>,
    <S::MdsArray as Index<usize>>::Output: Index<usize, Output = Zq>,
{
    pub fn new(input: Vec<Zq>, permutation: PoseidonPermutation<S>) -> Self {
        Self {
            input,
            output: vec![],
            permutation,
        }
    }

    // Absorb input into the sponge
    fn absorb(&mut self) -> Result<(), PoseidonError> {
        let mut remaining = self.input.clone();
        let mut position = 0;

        while !remaining.is_empty() {
            let available_space = S::RATE - position;
            let elements_to_absorb = remaining.len().min(available_space);
            let elements_to_process = &remaining[..elements_to_absorb];

            for (i, elem) in elements_to_process.iter().enumerate() {
                let target_index = S::CAPACITY + position + i;
                if target_index >= self.permutation.state.len() {
                    return Err(PoseidonError::OutOfBounds(target_index));
                }
                self.permutation.state[target_index] += *elem;
            }

            remaining = (&remaining[elements_to_absorb..]).to_vec();

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
    pub fn squeez(&mut self) -> Result<Vec<Zq>, SpongeError> {
        // Check if the requested number of elements exceeds the rate.
        if S::OUTPUT_LENGTH > S::RATE {
            return Err(SpongeError::OversizedSqueeze {
                requested: S::OUTPUT_LENGTH,
                rate: S::RATE,
            });
        }

        let mut output = vec![Zq::ZERO; S::OUTPUT_LENGTH];
        let mut remaining_elements = output.as_mut_slice();
        let mut squeeze_index = 0;
        // Loop until we return the expected amount
        loop {
            if squeeze_index + remaining_elements.len() <= S::RATE {
                remaining_elements.copy_from_slice(
                    &self.permutation.state[S::CAPACITY + squeeze_index
                        ..S::CAPACITY + remaining_elements.len() + squeeze_index],
                );
                return Ok(output);
            }
            // Reset squeeze index in case we need to squeeze again to fill the output
            let num_elements_squeezed = S::RATE - squeeze_index;
            remaining_elements[..num_elements_squeezed].copy_from_slice(
                &self.permutation.state[S::CAPACITY + squeeze_index
                    ..S::CAPACITY + num_elements_squeezed + squeeze_index],
            );

            remaining_elements = &mut remaining_elements[num_elements_squeezed..];
            self.permutation.permute();
            squeeze_index = 0;
        }
    }

    pub fn compute_hash(&mut self) -> &Vec<Zq> {
        if self.output.is_empty() {
            let _ = self.absorb();
            let _ = self.squeez();
        }
        &self.output
    }

    fn reset(&mut self) {
        self.input = vec![];
        self.output = vec![];
        self.permutation = PoseidonPermutation::new();
    }
}

// #[derive(Debug, thiserror::Error)]
// pub enum SpongeError {
//     #[error("Squeeze request of {requested} exceeds rate {rate}")]
//     OversizedSqueeze { requested: usize, rate: usize },
// }
// #[derive(Debug)]
// pub enum PoseidonError {
//     OutOfBounds(usize),
//     EmptyInput,
// }

// pub struct PoseidonHash<
//     const RATE: usize,
//     const CAPACITY: usize,
//     const WIDTH: usize,          // Width = Rate + CAPACiTY
//     const FULL_ROUNDS: usize, // Number of full rounds (S-box is applied to every element of the state)
//     const PARTIAL_ROUNDS: usize, // Number of partial rounds (S-box is applied to first element of the state)\
//     const ROUNDS: usize,         // Rounds = Full_rounds + Partial_rounds
//     const ALPHA: u64,
//     const OUTPUT_LENGTH: usize,
// > {
//     permutation: PoseidonPermutation<WIDTH, FULL_ROUNDS, PARTIAL_ROUNDS, ROUNDS, ALPHA>,
//     input: Vec<Zq>,
//     output: [Zq; OUTPUT_LENGTH],
// }

// impl<
//         const RATE: usize,
//         const CAPACITY: usize,
//         const WIDTH: usize,
//         const FULL_ROUNDS: usize,
//         const PARTIAL_ROUNDS: usize,
//         const ROUNDS: usize,
//         const ALPHA: u64,
//         const OUTPUT_LENGTH: usize,
//     >
//     PoseidonHash<RATE, CAPACITY, WIDTH, FULL_ROUNDS, PARTIAL_ROUNDS, ROUNDS, ALPHA, OUTPUT_LENGTH>
// {
//     /// Create a new Poseidon hash instance from the given input vector.
//     pub fn new(input: Vec<Zq>) -> Self {
//         let permutation_setup = PoseidonPermutation::<
//             WIDTH,
//             FULL_ROUNDS,
//             PARTIAL_ROUNDS,
//             ROUNDS,
//             ALPHA,
//         >::new(
//             PoseidonPermutation::<WIDTH, FULL_ROUNDS, PARTIAL_ROUNDS, ROUNDS, ALPHA>::generate_mds(
//             ),
//             PoseidonPermutation::<WIDTH, FULL_ROUNDS, PARTIAL_ROUNDS, ROUNDS, ALPHA>::generate_ark(
//             ),
//         );
//         Self {
//             permutation: permutation_setup,
//             input,
//             output: [Zq::ZERO; OUTPUT_LENGTH],
//         }
//     }

//     // Absorb input into the sponge
//     pub fn absorb(&mut self) -> Result<(), PoseidonError> {
//         let mut remaining = self.input.clone();
//         let mut position = 0;

//         while !remaining.is_empty() {
//             let available_space = RATE - position;
//             let elements_to_absorb = remaining.len().min(available_space);
//             let elements_to_process = &remaining[..elements_to_absorb];

//             for (i, elem) in elements_to_process.iter().enumerate() {
//                 let target_index = CAPACITY + position + i;
//                 if target_index >= self.permutation.state.len() {
//                     return Err(PoseidonError::OutOfBounds(target_index));
//                 }
//                 self.permutation.state[target_index] += *elem;
//             }

//             remaining = (&remaining[elements_to_absorb..]).to_vec();

//             if remaining.is_empty() {
//                 self.permutation.permute();
//                 return Ok(());
//             }
//             self.permutation.permute();
//             position = 0;
//         }

//         Ok(())
//     }

//     // Squeeze output from the sponge
//     pub fn squeeze(&mut self) -> Result<Vec<Zq>, SpongeError> {
//         // Check if the requested number of elements exceeds the rate.
//         if OUTPUT_LENGTH > RATE {
//             return Err(SpongeError::OversizedSqueeze {
//                 requested: OUTPUT_LENGTH,
//                 rate: RATE,
//             });
//         }

//         let mut output = vec![Zq::ZERO; OUTPUT_LENGTH];
//         let mut remaining_elements = output.as_mut_slice();
//         let mut squeeze_index = 0;
//         // Loop until we return the expected amount
//         loop {
//             if squeeze_index + remaining_elements.len() <= RATE {
//                 remaining_elements.copy_from_slice(
//                     &self.permutation.state[CAPACITY + squeeze_index
//                         ..CAPACITY + remaining_elements.len() + squeeze_index],
//                 );
//                 return Ok(output);
//             }
//             // Reset squeeze index in case we need to squeeze again to fill the output
//             let num_elements_squeezed = RATE - squeeze_index;
//             remaining_elements[..num_elements_squeezed].copy_from_slice(
//                 &self.permutation.state
//                     [CAPACITY + squeeze_index..CAPACITY + num_elements_squeezed + squeeze_index],
//             );

//             remaining_elements = &mut remaining_elements[num_elements_squeezed..];
//             self.permutation.permute();
//             squeeze_index = 0;
//         }
//     }

//     fn reset(&mut self) {
//         todo!()
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::zq::Zq;

//     #[test]
//     fn test_poseidon_hash_new_initialises_state_and_output() {
//         // Arrange
//         type Hash = PoseidonHash<2, 4, 6, 6, 35, 41, 3, 2>;
//         let input = vec![Zq::new(17), Zq::new(23)];

//         // Act
//         let hash = Hash::new(input);

//         // Assert
//         for cell in hash.output.iter() {
//             assert_eq!(*cell, Zq::ZERO);
//         }
//         assert_eq!(hash.permutation.ark.len(), 41); // ROUNDS
//         assert_eq!(hash.input.len(), 2);
//         assert_eq!(hash.output.len(), 2);
//         assert_eq!(hash.permutation.ark[0].len(), 6); // WIDTH
//         assert_eq!(hash.permutation.mds.len(), 6); // WIDTH × WIDTH matrix
//         assert_eq!(hash.permutation.mds[0].len(), 6);
//     }
// }
