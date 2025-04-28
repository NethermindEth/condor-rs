use crate::zq::Zq;
use blake2::Blake2b256;
use blake2::Digest;
use rand::distr::{Distribution, Uniform};
use rand::{random, rng, CryptoRng};

pub struct PoseidonPermutation<
    const WIDTH: usize,
    const ROUNDS: usize,
    const PARTIAL_ROUNDS: usize,
    const ALPHA: u64,
> {
    pub state: [Zq; WIDTH],
    pub ark: [[Zq; WIDTH]; ROUNDS],
    pub mds: [[Zq; WIDTH]; WIDTH],
}

impl<const WIDTH: usize, const ROUNDS: usize, const PARTIAL_ROUNDS: usize, const ALPHA: u64>
    PoseidonPermutation<WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>
{
    pub fn new() -> Self {
        Self {
            state: [Zq::ZERO; WIDTH],
            ark: Self::generate_round_constants(),
            mds: Self::generate_mds_matrix(rng()),
        }
    }

    pub fn new_with_ark_mds(ark: [[Zq; WIDTH]; ROUNDS], mds: [[Zq; WIDTH]; WIDTH]) -> Self {
        Self {
            state: [Zq::ZERO; WIDTH],
            ark,
            mds,
        }
    }

    // Generates ARK matrix
    fn generate_round_constants() -> [[Zq; WIDTH]; ROUNDS] {
        // Create a matrix to store constants for each round and state position
        let mut ark = [[Zq::new(0); WIDTH]; ROUNDS];

        // Iterate over the rows (rounds)
        for row in ark.iter_mut().take(ROUNDS) {
            // Iterate over the states in the current row using `iter_mut` for mutable access
            for state in row.iter_mut().take(WIDTH) {
                // 4 random bytes to use as input (seed)
                let seed: [u8; 4] = random();

                // Hasher based on Blake2b
                let mut hasher = Blake2b256::new();
                hasher.update(seed);

                // Get the hash output
                let hash_result: [u8; 32] = hasher.finalize().into();

                // Take the first 4 bytes of the hash
                let hash_bytes: [u8; 4] = hash_result[..4]
                    .try_into()
                    .expect("can cast a slice into an array");

                // Convert those 4 bytes into a u32 value
                let constant = u32::from_le_bytes(hash_bytes);

                // Store the constant in the matrix
                *state = Zq::new(constant);
            }
        }
        ark
    }

    /// Generates an MDS (Maximum Distance Separable) matrix using a Cauchy matrix
    fn generate_mds_matrix<R: CryptoRng>(mut rng: R) -> [[Zq; WIDTH]; WIDTH] {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::MAX).unwrap();

        let mut x_vals = Vec::new();
        let mut y_vals = Vec::new();

        // Sample unique x values
        while x_vals.len() < WIDTH {
            let candidate = uniform.sample(&mut rng);
            if !x_vals.contains(&candidate) {
                x_vals.push(candidate);
            }
        }

        // Sample unique y values with additional constraint
        while y_vals.len() < WIDTH {
            let candidate = uniform.sample(&mut rng);
            let valid =
                !y_vals.contains(&candidate) && !x_vals.iter().any(|x| *x + candidate == Zq::ZERO);
            if valid {
                y_vals.push(candidate);
            }
        }

        // Construct the Cauchy matrix
        let mut matrix = [[Zq::ZERO; WIDTH]; WIDTH];

        for (i, &x_val) in x_vals.iter().enumerate() {
            for (j, &y_val) in y_vals.iter().enumerate() {
                let denom = x_val + y_val;
                matrix[i][j] = denom.inv().expect("Inverse must exist here");
            }
        }

        matrix
    }

    fn apply_round_constants(&mut self, round_num: usize) {
        for (i, elem) in self.state.iter_mut().enumerate() {
            *elem += self.ark[round_num][i];
        }
    }

    fn apply_s_box(&mut self, is_full_round: bool) {
        if is_full_round {
            for elem in &mut self.state {
                // Apply transformation to each element of state
                *elem = elem.pow(ALPHA);
            }
        } else {
            // Apply transformation only to the first element
            self.state[0] = self.state[0].pow(ALPHA);
        }
    }

    fn apply_mds(&mut self) {
        let mut old_state = self.state.clone();
        // Matrix multiplication with the state for diffusion
        for (i, _cur) in self.state.iter().enumerate() {
            let mut sum = Zq::ZERO;
            for (j, elem) in self.state.iter().enumerate() {
                // matrix multiplication
                sum += *elem * self.mds[i][j];
            }
            old_state[i] = sum;
        }
        self.state = old_state;
    }

    pub fn permute(&mut self) {
        let full_rounds_half = (ROUNDS - PARTIAL_ROUNDS) / 2;
        // full round ark/S-box/mds
        for i in 0..full_rounds_half {
            self.apply_round_constants(i);
            self.apply_s_box(true);
            self.apply_mds();
        }
        // partial round ark/S-box/mds
        for i in full_rounds_half..(full_rounds_half + PARTIAL_ROUNDS) {
            self.apply_round_constants(i);
            self.apply_s_box(false); // Apply the S-box to the first element only
            self.apply_mds();
        }
        // full round ark/S-box/mds
        for i in (full_rounds_half + PARTIAL_ROUNDS)..(ROUNDS) {
            self.apply_round_constants(i);
            self.apply_s_box(true);
            self.apply_mds();
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::poseidon::permutation::PoseidonPermutation;
    use crate::zq::Zq;

    type Permutation1 = PoseidonPermutation<6, 41, 6, 3>;
    type Permutation2 = PoseidonPermutation<6, 6, 6, 3>;
    type Permutation3 = PoseidonPermutation<2, 2, 0, 1>;
    type Permutation4 = PoseidonPermutation<2, 3, 1, 2>;

    #[test]
    fn test_correct_permutation_config() {
        let permutation = Permutation1::new();
        assert_eq!(permutation.state.len(), 6);
    }

    #[test]
    fn test_s_box_partial_round_alpha_3() {
        // Arrange
        let mut permutation = Permutation1::new();
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
        let mut permutation = Permutation1::new();
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
        let mut permutation = Permutation2::new();
        permutation.ark = [
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
        let mut permutation = Permutation1::new();
        permutation.mds = [
            [Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE],
            [Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE],
            [Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE],
            [Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE],
            [Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE],
            [Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE, Zq::ONE],
        ];
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
        let mut permutation = Permutation3::new();
        permutation.ark = [[Zq::new(5), Zq::new(7)], [Zq::new(11), Zq::new(13)]];
        permutation.mds = [[Zq::new(1), Zq::new(0)], [Zq::new(0), Zq::new(1)]];

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
        let mut permutation = Permutation4::new();
        permutation.ark = [[Zq::ZERO; 2]; 3];
        permutation.mds = [[Zq::new(1), Zq::new(0)], [Zq::new(0), Zq::new(1)]];
        permutation.state[0] = Zq::new(2);
        permutation.state[1] = Zq::new(3);

        // Act
        permutation.permute();

        // Assert
        assert_eq!(permutation.state[0], Zq::new(256)); // 2^(2→4→16→256)
        assert_eq!(permutation.state[1], Zq::new(81)); // 3^(2→9→9→81)
    }
}
