// use crate::zq::Zq;
// use blake2::Blake2b256;
// use blake2::Digest;
// use rand::distr::{Distribution, Uniform};
// use rand::{random, CryptoRng};


// use crate::poseidon::sponge::Pose;
// use crate::poseidon::sponge::{PoseidonError, SpongeError};
// use crate::zq::Zq;
// use blake2::Blake2b256;
// use blake2::Digest;
// use rand::distr::{Distribution, Uniform};
// use rand::{random, CryptoRng};
// use std::convert::TryInto;

// pub trait Transcript {
//     fn new(rng: &mut impl CryptoRng) -> Self;
//     fn absorb(&mut self, value: Zq) -> Result<(), PoseidonError>;
//     fn get_challenge(&mut self) -> Result<Zq, SpongeError>;
// }


// pub struct PoseidonTranscript {
//     sponge: Pose,
// }

// impl Transcript for PoseidonTranscript {
//     fn new(rng: &mut impl CryptoRng) -> Self {
//         let sponge = crate::poseidon::PoseidonSponge::new(
//             8,
//             8,
//             5,
//             8,
//             8,
//             cauchy_mds_matrix(rng, 16),
//             generate_secure_round_constants(13, 16),
//         );

//         Self { sponge }
//     }

//     fn absorb(&mut self, value: Zq) -> Result<(), PoseidonError> {
//         self.sponge.absorb(&[value])?; // ? Propagates the error if it occurs
//         Ok(())
//     }

//     fn get_challenge(&mut self) -> Result<Zq, SpongeError> {
//         let squeezed_elements = self.sponge.squeeze(1);

//         match squeezed_elements {
//             Ok(vec) => Ok(vec[0] + Zq::ONE), // // Example: return first challenge + 1
//             Err(err) => Err(err),            // Propagate the SpongeError
//         }
//     }
// }
