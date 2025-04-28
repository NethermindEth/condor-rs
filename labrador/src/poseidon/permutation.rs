use crate::poseidon::constants::PermutationConstants;
use crate::zq::Zq;
use std::marker::PhantomData;
use std::ops::Index;

#[derive(Default)]
pub struct PoseidonPermutation<P: PermutationConstants> {
    pub state: Vec<Zq>,
    _marker: PhantomData<P>,
}

impl<P> PoseidonPermutation<P>
where
    P: PermutationConstants,
    P::ArkArray: Index<usize>,
    P::MdsArray: Index<usize>,
    <P::ArkArray as Index<usize>>::Output: Index<usize, Output = Zq>,
    <P::MdsArray as Index<usize>>::Output: Index<usize, Output = Zq>,
{
    pub fn new() -> Self {
        Self {
            state: vec![Zq::ZERO; P::WIDTH],
            _marker: PhantomData,
        }
    }

    pub fn apply_round_constants(&mut self, round_num: usize) {
        for (i, elem) in self.state.iter_mut().enumerate() {
            *elem += P::ARK[round_num][i];
        }
    }

    pub fn apply_s_box(&mut self, is_full_round: bool) {
        if is_full_round {
            for elem in &mut self.state {
                // Apply transformation to each element of state
                *elem = elem.pow(P::ALPHA);
            }
        } else {
            // Apply transformation only to the first element
            self.state[0] = self.state[0].pow(P::ALPHA);
        }
    }

    pub fn apply_mds(&mut self) {
        let mut old_state = self.state.clone();
        // Matrix multiplication with the state for diffusion
        for (i, _cur) in self.state.iter().enumerate() {
            let mut sum = Zq::ZERO;
            for (j, elem) in self.state.iter().enumerate() {
                // matrix multiplication
                sum += *elem * P::MDS[i][j];
            }
            old_state[i] = sum;
        }
        self.state = old_state;
    }

    pub fn permute(&mut self) {
        let full_rounds_half = P::FULL_ROUNDS / 2;
        // full round ark/S-box/mds
        for i in 0..full_rounds_half {
            self.apply_round_constants(i);
            self.apply_s_box(true);
            self.apply_mds();
        }
        // partial round ark/S-box/mds
        for i in full_rounds_half..(full_rounds_half + P::PARTIAL_ROUNDS) {
            self.apply_round_constants(i);
            self.apply_s_box(false); // Apply the S-box to the first element only
            self.apply_mds();
        }
        // full round ark/S-box/mds
        for i in (full_rounds_half + P::PARTIAL_ROUNDS)..(P::PARTIAL_ROUNDS + P::FULL_ROUNDS) {
            self.apply_round_constants(i);
            self.apply_s_box(true);
            self.apply_mds();
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::poseidon::constants::PermutationConstants;
    use crate::poseidon::permutation::PoseidonPermutation;
    use crate::zq::Zq;

    struct DummyPermutation;
    impl PermutationConstants for DummyPermutation {
        const WIDTH: usize = 6;
        const PARTIAL_ROUNDS: usize = 6;
        const FULL_ROUNDS: usize = 35;
        const ALPHA: u64 = 3;

        type ArkArray = [[Zq; Self::WIDTH]; Self::PARTIAL_ROUNDS + Self::FULL_ROUNDS]; // Concrete sizes required here
        type MdsArray = [[Zq; Self::WIDTH]; Self::WIDTH]; // Concrete sizes required here

        const ARK: Self::ArkArray =
            [[Zq::new(0); Self::WIDTH]; Self::PARTIAL_ROUNDS + Self::FULL_ROUNDS];
        const MDS: Self::MdsArray = [[Zq::new(1); Self::WIDTH]; Self::WIDTH];

        fn ark() -> Self::ArkArray {
            Self::ARK
        }

        fn mds() -> Self::MdsArray {
            Self::MDS
        }
    }

    #[test]
    fn test_correct_permutation_config() {
        let permutation = PoseidonPermutation::<DummyPermutation>::new();
        assert_eq!(permutation.state.len(), 6);
    }

    #[test]
    fn test_s_box_partial_round_alpha_3() {
        // Arrange
        let mut permutation = PoseidonPermutation::<DummyPermutation>::new();
        for (ctr, elem) in permutation.state.iter_mut().enumerate() {
            *elem = Zq::new(ctr as u32 + 10);
        }

        // Act
        permutation.apply_s_box(false);

        // Assert
        assert_eq!(permutation.state[0], Zq::new(10).pow(3));
        for i in 1..(permutation.state.len()) {
            assert_eq!(permutation.state[i], Zq::new(i as u32 + 10));
        }
    }

    #[test]
    fn test_s_box_full_round_alpha_3() {
        // Arrange
        let mut permutation = PoseidonPermutation::<DummyPermutation>::new();
        for (ctr, elem) in permutation.state.iter_mut().enumerate() {
            *elem = Zq::new(ctr as u32 + 10);
        }

        // Act
        permutation.apply_s_box(true);

        // Assert
        for (index, elem) in permutation.state.iter().enumerate() {
            assert_eq!(elem, &Zq::new(index as u32 + 10).pow(3));
        }
    }

    #[test]
    fn test_apply_ark() {
        // Arrange
        struct ARKPermutation;
        impl PermutationConstants for ARKPermutation {
            const WIDTH: usize = 6;
            const PARTIAL_ROUNDS: usize = 6;
            const FULL_ROUNDS: usize = 0;
            const ALPHA: u64 = 3;

            type ArkArray = [[Zq; Self::WIDTH]; Self::PARTIAL_ROUNDS + Self::FULL_ROUNDS]; // Concrete sizes required here
            type MdsArray = [[Zq; Self::WIDTH]; Self::WIDTH]; // Concrete sizes required here

            const ARK: Self::ArkArray = [
                [
                    Zq::new(0),
                    Zq::new(1),
                    Zq::new(2),
                    Zq::new(3),
                    Zq::new(4),
                    Zq::new(5),
                ],
                [Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO],
                [Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO],
                [Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO],
                [Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO],
                [Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO],
            ];

            const MDS: Self::MdsArray = [[Zq::ZERO; Self::WIDTH]; Self::WIDTH];
        }

        let mut permutation = PoseidonPermutation::<ARKPermutation>::new();

        // Act
        permutation.apply_round_constants(0);

        // Assert
        for (idx, elem) in permutation.state.iter().enumerate() {
            assert_eq!(elem, &Zq::new(idx as u32));
        }
    }

    #[test]
    fn test_apply_mds_all_one_matrix() {
        // Arrange
        let mut permutation = PoseidonPermutation::<DummyPermutation>::new();
        for (i, elem) in permutation.state.iter_mut().enumerate() {
            *elem = Zq::new((i + 1) as u32); // state = [1,2,3,4,5,6]
        }

        // Act
        permutation.apply_mds();

        // Assert
        let expected = Zq::new(21);
        for elem in &permutation.state {
            assert_eq!(elem, &expected);
        }
    }

    /// permute: ALPHA = 1 (identity S-box) and identity MDS.
    /// With two full rounds and *no* partial rounds the permutation
    /// should simply add the round constants twice.
    #[test]
    fn test_permute_with_identity_mds_s_box() {
        // Arrange
        struct SBoxPermutation;
        impl PermutationConstants for SBoxPermutation {
            const WIDTH: usize = 2;
            const PARTIAL_ROUNDS: usize = 0;
            const FULL_ROUNDS: usize = 2;
            const ALPHA: u64 = 1;

            type ArkArray = [[Zq; Self::WIDTH]; Self::PARTIAL_ROUNDS + Self::FULL_ROUNDS]; // Concrete sizes required here
            type MdsArray = [[Zq; Self::WIDTH]; Self::WIDTH]; // Concrete sizes required here

            const ARK: Self::ArkArray = [[Zq::new(5), Zq::new(7)], [Zq::new(11), Zq::new(13)]];
            const MDS: Self::MdsArray = [[Zq::new(1), Zq::new(0)], [Zq::new(0), Zq::new(1)]];
        }
        let mut permutation = PoseidonPermutation::<SBoxPermutation>::new();

        // Act
        permutation.permute();

        // Assert
        assert_eq!(permutation.state[0], Zq::new(16)); // 5+11
        assert_eq!(permutation.state[1], Zq::new(20)); // 7+13
    }

    /// permute: ALPHA = 2 (square), **one partial round**.
    /// Identity MDS and zero ARK so the result depends *only*
    /// on the S-box schedule:
    ///
    ///   initial:     (2,3)
    ///   round-0 F:   ( 4,  9)
    ///   round-1 P:   (16,  9)
    ///   round-2 F:   (256, 81)
    #[test]
    fn test_permute_partial_round_sbox() {
        // Arrange
        struct SBoxPermutation;
        impl PermutationConstants for SBoxPermutation {
            const WIDTH: usize = 2;
            const PARTIAL_ROUNDS: usize = 1;
            const FULL_ROUNDS: usize = 2;
            const ALPHA: u64 = 2;

            type ArkArray = [[Zq; Self::WIDTH]; Self::PARTIAL_ROUNDS + Self::FULL_ROUNDS]; // Concrete sizes required here
            type MdsArray = [[Zq; Self::WIDTH]; Self::WIDTH]; // Concrete sizes required here

            const ARK: Self::ArkArray = [[Zq::ZERO; 2]; 3];
            const MDS: Self::MdsArray = [[Zq::new(1), Zq::new(0)], [Zq::new(0), Zq::new(1)]];
        }
        let mut permutation = PoseidonPermutation::<SBoxPermutation>::new();
        permutation.state[0] = Zq::new(2);
        permutation.state[1] = Zq::new(3);

        // Act
        permutation.permute();

        // Assert
        assert_eq!(permutation.state[0], Zq::new(256)); // 2^(2→4→16→256)
        assert_eq!(permutation.state[1], Zq::new(81)); // 3^(2→9→9→81)
    }
}

// pub struct PoseidonPermutation<
//     const WIDTH: usize,          // Width = Rate + Capacity
//     const FULL_ROUNDS: usize, // Number of full rounds (S-box is applied to every element of the state)
//     const PARTIAL_ROUNDS: usize, // Number of partial rounds (S-box is applied to first element of the state)\
//     const ROUNDS: usize,         // Rounds = Full_rounds + Partial_rounds
//     const ALPHA: u64,            // Exponent for the S-box
// > {
//     pub state: [Zq; WIDTH],         // Current state of the sponge
//     pub mds: [[Zq; WIDTH]; WIDTH],  // MDS matrix [[Zq; RATE + CAPACITY]; RATE + CAPACITY]
//     pub ark: [[Zq; WIDTH]; ROUNDS], // Round constants [[Zq; RATE + CAPACITY]; FULL_ROUNDS + PARTIAL_ROUNDS]
// }

// impl<
//         const WIDTH: usize,
//         const FULL_ROUNDS: usize,
//         const PARTIAL_ROUNDS: usize,
//         const ROUNDS: usize,
//         const ALPHA: u64,
//     > PoseidonPermutation<WIDTH, FULL_ROUNDS, PARTIAL_ROUNDS, ROUNDS, ALPHA>
// {
//     fn generate_state() -> [Zq; WIDTH] {
//         [Zq::ZERO; WIDTH]
//     }

//     pub fn generate_mds() -> [[Zq; WIDTH]; WIDTH] {
//         [[Zq::default(); WIDTH]; WIDTH]
//     }

//     pub fn generate_ark() -> [[Zq; WIDTH]; ROUNDS] {
//         [[Zq::default(); WIDTH]; ROUNDS]
//     }

//     // Initialize Poseidon sponge parameters
//     pub fn new(mds: [[Zq; WIDTH]; WIDTH], ark: [[Zq; WIDTH]; ROUNDS]) -> Self {
//         Self {
//             state: Self::generate_state(),
//             mds,
//             ark,
//         }
//     }

//     // Apply the round constants (ARK)
//     fn apply_ark(&mut self, round_num: usize) {
//         for (i, elem) in self.state.iter_mut().enumerate() {
//             // adds diffusion (non-linearity)
//             *elem += self.ark[round_num][i];
//         }
//     }

//     // Apply the S-box (x -> x^alpha)
//     fn apply_s_box(&mut self, is_full_round: bool) {
//         if is_full_round {
//             for elem in &mut self.state {
//                 // Apply transformation to each element of state
//                 *elem = elem.pow(ALPHA);
//             }
//         } else {
//             // Apply transformation only to the first element
//             self.state[0] = self.state[0].pow(ALPHA);
//         }
//     }

//     // Apply the MDS (Maximum Distance Separable) matrix
//     fn apply_mds(&mut self) {
//         let mut new_state = self.state.clone();

//         // Matrix multiplication with the state for diffusion
//         for (i, _cur) in self.state.iter().enumerate() {
//             let mut sum = Zq::ZERO;
//             for (j, elem) in self.state.iter().enumerate() {
//                 // matrix multiplication
//                 sum += *elem * self.mds[i][j];
//             }
//             new_state[i] = sum;
//         }

//         self.state = new_state;
//     }

//     // Apply a permutation
//     pub fn permute(&mut self) {
//         let full_rounds_half = FULL_ROUNDS / 2;

//         // full round ark/S-box/mds
//         for i in 0..full_rounds_half {
//             self.apply_ark(i);
//             self.apply_s_box(true);
//             self.apply_mds();
//         }
//         // partial round ark/S-box/mds
//         for i in full_rounds_half..(full_rounds_half + PARTIAL_ROUNDS) {
//             self.apply_ark(i);
//             self.apply_s_box(false); // Apply the S-box to the first element only
//             self.apply_mds();
//         }
//         // full round ark/S-box/mds
//         for i in (full_rounds_half + PARTIAL_ROUNDS)..(PARTIAL_ROUNDS + FULL_ROUNDS) {
//             self.apply_ark(i);
//             self.apply_s_box(true);
//             self.apply_mds();
//         }
//     }
// }

// #[cfg(test)]
// mod tests {
//     use crate::poseidon::permutation::PoseidonPermutation;
//     use crate::zq::Zq;

//     #[test]
//     fn test_correct_array_length() {
//         type SpongeType = PoseidonPermutation<6, 6, 35, 41, 3>;
//         let sponge = SpongeType::new([[Zq::new(1); 6]; 6], [[Zq::new(1); 6]; 41]);
//         assert_eq!(sponge.state.len(), 6);
//         assert_eq!(sponge.ark.len(), 41);
//         assert_eq!(sponge.ark[0].len(), 6);
//         assert_eq!(sponge.mds.len(), 6);
//         assert_eq!(sponge.mds[0].len(), 6);
//     }

//     #[test]
//     fn test_s_box_partial_round_alpha_3() {
//         // Arrange
//         type SpongeType = PoseidonPermutation<6, 6, 35, 41, 3>;
//         let mut sponge = SpongeType::new([[Zq::default(); 6]; 6], [[Zq::new(1); 6]; 41]);
//         for (ctr, elem) in sponge.state.iter_mut().enumerate() {
//             *elem = Zq::new(ctr as u32 + 10);
//         }

//         // Act
//         sponge.apply_s_box(false);

//         // Assert
//         assert_eq!(sponge.state[0], Zq::new(10).pow(3));
//         for i in 1..(sponge.state.len()) {
//             assert_eq!(sponge.state[i], Zq::new(i as u32 + 10));
//         }
//     }

//     #[test]
//     fn test_s_box_full_round_alpha_3() {
//         // Arrange
//         type SpongeType = PoseidonPermutation<6, 6, 35, 41, 3>;
//         let mut sponge = SpongeType::new([[Zq::default(); 6]; 6], [[Zq::new(1); 6]; 41]);
//         for (ctr, elem) in sponge.state.iter_mut().enumerate() {
//             *elem = Zq::new(ctr as u32 + 10);
//         }

//         // Act
//         sponge.apply_s_box(true);

//         // Assert
//         for (index, elem) in sponge.state.iter().enumerate() {
//             assert_eq!(elem, &Zq::new(index as u32 + 10).pow(3));
//         }
//     }

//     #[test]
//     fn test_apply_ark() {
//         // Arrange
//         type SpongeType = PoseidonPermutation<6, 6, 35, 41, 3>;
//         let mut ark = [[Zq::ZERO; 6]; 41];
//         for (idx, rc) in ark[0].iter_mut().enumerate() {
//             *rc = Zq::new((idx as u32) + 1);
//         }
//         let mut permutation = SpongeType::new([[Zq::default(); 6]; 6], ark);

//         // Act
//         permutation.apply_ark(0);

//         // Assert
//         for (idx, elem) in permutation.state.iter().enumerate() {
//             assert_eq!(elem, &Zq::new((idx as u32) + 1));
//         }
//     }

//     #[test]
//     fn test_apply_mds_all_one_matrix() {
//         // Arrange
//         type SpongeType = PoseidonPermutation<6, 6, 35, 41, 3>;
//         let mds = [[Zq::new(1); 6]; 6];
//         let ark = [[Zq::ZERO; 6]; 41];
//         let mut sponge = SpongeType::new(mds, ark);
//         for (i, elem) in sponge.state.iter_mut().enumerate() {
//             *elem = Zq::new((i + 1) as u32); // state = [1,2,3,4,5,6]
//         }

//         // Act
//         sponge.apply_mds();

//         // Assert
//         let expected = Zq::new(21);
//         for elem in &sponge.state {
//             assert_eq!(elem, &expected);
//         }
//     }

//     /// permute: ALPHA = 1 (identity S-box) and identity MDS.
//     /// With two full rounds and *no* partial rounds the permutation
//     /// should simply add the round constants twice.
//     #[test]
//     fn test_permute_with_identity_mds_s_box() {
//         // Arrange
//         type Sponge = PoseidonPermutation<2, 2, 0, 2, 1>; // WIDTH=2, F=2, P=0
//                                                           // Identity MDS
//         let mds = [[Zq::new(1), Zq::new(0)], [Zq::new(0), Zq::new(1)]];
//         // Round-0 = [5,7] ; Round-1 = [11,13]
//         let mut ark = [[Zq::ZERO; 2]; 2];
//         ark[0] = [Zq::new(5), Zq::new(7)];
//         ark[1] = [Zq::new(11), Zq::new(13)];
//         let mut sponge = Sponge::new(mds, ark);

//         // Act
//         sponge.permute();

//         // Assert
//         assert_eq!(sponge.state[0], Zq::new(16)); // 5+11
//         assert_eq!(sponge.state[1], Zq::new(20)); // 7+13
//     }

//     /// permute: ALPHA = 2 (square), **one partial round**.
//     /// Identity MDS and zero ARK so the result depends *only*
//     /// on the S-box schedule:
//     ///
//     ///   initial:     (2,3)
//     ///   round-0 F:   ( 4,  9)
//     ///   round-1 P:   (16,  9)
//     ///   round-2 F:   (256, 81)
//     #[test]
//     fn test_permute_partial_round_sbox() {
//         // Arrange
//         type Sponge = PoseidonPermutation<2, 2, 1, 3, 2>;
//         let mds = [[Zq::new(1), Zq::new(0)], [Zq::new(0), Zq::new(1)]];
//         let ark = [[Zq::ZERO; 2]; 3]; // no round constants
//         let mut sponge = Sponge::new(mds, ark);
//         sponge.state[0] = Zq::new(2);
//         sponge.state[1] = Zq::new(3);

//         // Act
//         sponge.permute();

//         // Assert
//         assert_eq!(sponge.state[0], Zq::new(256)); // 2^(2→4→16→256)
//         assert_eq!(sponge.state[1], Zq::new(81)); // 3^(2→9→9→81)
//     }
// }
