pub trait PermutationConstants {
    const WIDTH: usize;
    const PARTIAL_ROUNDS: usize;
    const FULL_ROUNDS: usize;
    const ALPHA: u64;

    type ArkArray;
    type MdsArray;
    const ARK: Self::ArkArray;
    const MDS: Self::MdsArray;

    fn ark() -> Self::ArkArray {
        Self::ARK
    }

    fn mds() -> Self::MdsArray {
        Self::MDS
    }
}

// #[test]
// fn test_correct_array_length() {
//     type SpongeType = PoseidonPermutation<6, 6, 35, 41, 3>;
//     let sponge = SpongeType::new([[Zq::new(1); 6]; 6], [[Zq::new(1); 6]; 41]);
//     assert_eq!(sponge.state.len(), 6);
//     assert_eq!(sponge.ark.len(), 41);
//     assert_eq!(sponge.ark[0].len(), 6);
//     assert_eq!(sponge.mds.len(), 6);
//     assert_eq!(sponge.mds[0].len(), 6);
// }

pub trait PoseidonHash {
    const RATE: usize;
    const CAPACITY: usize;
    const OUTPUT_LENGTH: usize;

    fn new();
    fn absorb();
    fn squeez();
    fn hash();
    fn reset();
}

// // Implementation example
// struct Poseidon3_8_30;

// impl PermutationSetup for Poseidon3_8_30 {
//     const WIDTH: usize = 3;
//     const PARTIAL_ROUNDS: usize = 30;
//     const FULL_ROUNDS: usize = 8;
//     const ROUNDS: usize = 38;  // 30 + 8
//     const ALPHA: u64 = 5;

//     type ArkArray = [[Zq; Self::WIDTH]; Self::ROUNDS];  // Now using concrete sizes
//     type MdsArray = [[Zq; 3]; 3];   // Now using concrete sizes

//     fn ark() -> Self::ArkArray {
//         // Implementation
//         [[Zq(0); 3]; 38]  // Placeholder
//     }

//     fn mds() -> Self::MdsArray {
//         // Implementation
//         [[Zq(0); 3]; 3]  // Placeholder
//     }
// }
