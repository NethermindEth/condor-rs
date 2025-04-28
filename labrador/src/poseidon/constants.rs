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

pub trait SpongeConstants: PermutationConstants {
    const RATE: usize;
    const CAPACITY: usize;
    const OUTPUT_LENGTH: usize;
}
