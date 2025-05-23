pub mod sponges;
use crate::{
    core::jl::Projection,
    ring::{rq::Rq, rq_matrix::RqMatrix, rq_vector::RqVector, zq::Zq},
};
pub use sponges::Sponge;

pub struct LabradorTranscript<S: Sponge> {
    sponge: S,
    pub u1: RqVector,
    pub vector_p: Vec<Zq>,
    pub b_ct_aggr: RqVector,
    pub u2: RqVector,
    pub z: RqVector,
    pub t: RqMatrix,
    pub g: RqMatrix,
    pub h: RqMatrix,
}

impl<S: Sponge> LabradorTranscript<S> {
    pub fn new(sponge: S) -> Self {
        Self {
            sponge,
            u1: RqVector::new(Vec::new()),
            vector_p: Vec::new(),
            b_ct_aggr: RqVector::new(Vec::new()),
            u2: RqVector::new(Vec::new()),
            z: RqVector::new(Vec::new()),
            t: RqMatrix::new(vec![RqVector::new(Vec::new())], false),
            g: RqMatrix::new(vec![RqVector::new(Vec::new())], true),
            h: RqMatrix::new(vec![RqVector::new(Vec::new())], true),
        }
    }

    pub fn set_u1(&mut self, u1: RqVector) {
        self.absorb_u1(&u1);
        self.u1 = u1;
    }

    pub fn absorb_u1(&mut self, u1: &RqVector) {
        self.sponge.absorb_rq(u1.get_elements());
    }

    pub fn set_vector_p(&mut self, p: Vec<Zq>) {
        self.absorb_vector_p(&p);
        self.vector_p = p;
    }

    pub fn absorb_vector_p(&mut self, p: &[Zq]) {
        self.sponge.absorb_zq(p);
    }

    pub fn set_vector_b_ct_aggr(&mut self, input: RqVector) {
        self.absorb_vector_b_ct_aggr(&input);
        self.b_ct_aggr = input;
    }

    pub fn absorb_vector_b_ct_aggr(&mut self, input: &RqVector) {
        self.sponge.absorb_rq(input.get_elements());
    }

    pub fn set_u2(&mut self, u2: RqVector) {
        self.absorb_u2(&u2);
        self.u2 = u2;
    }

    pub fn absorb_u2(&mut self, u2: &RqVector) {
        self.sponge.absorb_rq(u2.get_elements());
    }

    pub fn set_recursive_part(&mut self, z: RqVector, t: RqMatrix, g: RqMatrix, h: RqMatrix) {
        self.z = z;
        self.t = t;
        self.g = g;
        self.h = h;
    }

    pub fn generate_projections(
        &mut self,
        security_parameter: usize,
        rank: usize,
        multiplicity: usize,
    ) -> Projection {
        // r vectors, each of length 256 * nD
        let row_size = 2 * security_parameter;
        let col_size = rank * Rq::DEGREE;

        let matrices = (0..multiplicity)
            .map(|_| {
                let linear_projection_randomness = self.sponge.squeeze_zq(row_size * col_size);
                linear_projection_randomness
                    .chunks_exact(col_size)
                    .map(|chunk| {
                        let coeffs = chunk
                            .iter()
                            .map(|elem| {
                                if elem.get_value() < 2_u32.pow(30) {
                                    Zq::MAX
                                } else if elem.get_value() < 2_u32.pow(31) {
                                    Zq::ONE
                                } else {
                                    Zq::ZERO
                                }
                            })
                            .collect();
                        RqVector::new_from_zq_vector(coeffs)
                    })
                    .collect()
            })
            .collect();
        Projection::new(matrices, security_parameter)
    }

    pub fn generate_vector_psi(
        &mut self,
        number_of_vectors: usize,
        vector_length: usize,
    ) -> Vec<Vec<Zq>> {
        let elements = self.sponge.squeeze_zq(number_of_vectors * vector_length);
        elements
            .chunks_exact(vector_length)
            .map(|chunk| chunk.into())
            .collect()
    }

    pub fn generate_vector_omega(
        &mut self,
        number_of_vectors: usize,
        security_parameter: usize,
    ) -> Vec<Vec<Zq>> {
        let elements = self
            .sponge
            .squeeze_zq(number_of_vectors * 2 * security_parameter);
        elements
            .chunks_exact(2 * security_parameter)
            .map(|chunk| chunk.into())
            .collect()
    }

    pub fn generate_rq_vector(&mut self, vector_length: usize) -> RqVector {
        RqVector::new(self.sponge.squeeze_rq(vector_length))
    }

    pub fn generate_challenges(&mut self, op_norm: f64, multiplicity: usize) -> RqVector {
        let mut result = Vec::new();
        loop {
            let candidate = self.sample_challenge();
            if candidate.operator_norm() < op_norm {
                result.push(candidate);
                if result.len() >= multiplicity {
                    return RqVector::new(result);
                }
            }
        }
    }

    #[allow(clippy::as_conversions)]
    fn sample_challenge(&mut self) -> Rq {
        let mut coeffs = [Zq::ZERO; Rq::DEGREE];
        let random_bits = self.sponge.squeeze_bits(10 + 31);
        // Add 31 coefficients, each either +1 or -1.
        for (item, random_bit) in coeffs.iter_mut().zip(&random_bits).take(31) {
            *item = if *random_bit { Zq::new(1) } else { -Zq::new(1) };
        }
        // Add 10 coefficients, each either +2 or -2.
        for (item, random_bit) in coeffs.iter_mut().zip(random_bits).skip(31).take(10) {
            *item = if random_bit { Zq::new(2) } else { -Zq::new(2) };
        }

        // ------------------ 3. Shuffle using Fisher–Yates approach ---------
        for i in (1..Rq::DEGREE).rev() {
            let random_bytes = self.sponge.squeeze_bytes(4);
            let random_index = u32::from_le_bytes(random_bytes.try_into().unwrap());
            let random_index = (random_index as usize) % (i + 1); // uniform in 0..i
            coeffs.swap(i, random_index);
        }
        Rq::new(coeffs)
    }
}

#[cfg(test)]
mod tests_generate_pi {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_projection_matrix_has_correct_size() {
        let (security_parameter, rank, multiplicity) = (128, 20, 9);
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let projections = transcript.generate_projections(security_parameter, rank, multiplicity);
        assert_eq!(projections.get_projection_matrices().len(), multiplicity); // number_of_project_matrices
        assert_eq!(
            projections.get_projection_matrices()[0].get_row_len(),
            2 * security_parameter
        );
        assert_eq!(projections.get_projection_matrices()[0].get_col_len(), rank);
    }

    // Test the distribution of values in the random matrix
    #[test]
    #[allow(clippy::as_conversions)]
    fn test_projection_matrix_is_random() {
        let (security_parameter, rank, multiplicity) = (128, 1000, 1);
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let projections = transcript.generate_projections(security_parameter, rank, multiplicity);

        for projection_matrix in projections.get_projection_matrices() {
            let mut counts = [0.0, 0.0, 0.0]; // -1, 0, 1
            for row in projection_matrix.get_elements() {
                for cell in row.concatenate_coefficients() {
                    match cell {
                        Zq::ZERO => counts[1] += 1.0,
                        Zq::ONE => counts[2] += 1.0,
                        Zq::MAX => counts[0] += 1.0,
                        _ => panic!("Should not occur"),
                    }
                }
            }
            // Number of elements in the matrix as f64 (256x4x1000)
            #[allow(clippy::as_conversions)]
            let total: f64 = (256 * Rq::DEGREE * rank) as f64;
            println!("this is the total amount of elements{}", total);
            let expected = [0.25, 0.5, 0.25];
            for i in 0..3 {
                let actual = counts[i] / total;
                println!("This is the actual value {}", actual);
                assert!(
                    //Since its a statistical test some small error tolerance is allowed
                    (actual - expected[i]).abs() < 0.005,
                    "Values are not within expected proportions"
                );
            }
        }
    }
}

#[cfg(test)]
mod test_generate_psi {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_generate_vector_psi_has_correct_size() {
        let (k_range, l_range) = (20, 12);
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let result = transcript.generate_vector_psi(k_range, l_range);
        assert_eq!(result.len(), k_range); // number_of_project_matrices
        assert_eq!(result[0].len(), l_range);
    }

    // TODO: Testing randomness of the distribution is needed.
}

#[cfg(test)]
mod test_generate_omega {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_generate_vector_omega_has_correct_size() {
        let (security_parameter, k_range) = (128, 20);
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let result = transcript.generate_vector_omega(k_range, security_parameter);
        assert_eq!(result.len(), k_range); // number_of_project_matrices
        assert_eq!(result[0].len(), 256);
    }

    // TODO: Testing randomness of the distribution is needed.
}

#[cfg(test)]
mod test_generate_rq {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_generate_rq_has_correct_size() {
        let k_range = 20;
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let result = transcript.generate_rq_vector(k_range);
        assert_eq!(result.get_elements().len(), k_range); // number_of_project_matrices
    }

    // TODO: Testing randomness of the distribution is needed.
}

#[cfg(test)]
mod test_generate_challenges {
    use super::*;
    use crate::transcript::sponges::shake::ShakeSponge;

    #[test]
    fn test_generate_challenges_has_correct_size() {
        let multiplicity = 9;
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let result = transcript.generate_challenges(15.0, multiplicity);
        assert_eq!(result.get_elements().len(), multiplicity); // number_of_project_matrices
    }

    /// Test Challenge Set:
    /// const 71 and const 15 are from Challenge Space, paper page 6.
    /// l2 norm <= 71
    /// operator norm <= 15
    #[test]
    fn test_challenge_set() {
        let op_norm = 15.0;
        let multiplicity = 9;
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());

        let challenge_set = transcript.generate_challenges(15.0, multiplicity);

        for i in 0..multiplicity {
            // l2 norm 71 is from paper page 6, Challenge Space.
            assert_eq!(
                challenge_set[i].inner_product(&challenge_set[i]),
                Zq::new(71)
            );
            assert!(challenge_set[i].operator_norm() <= op_norm);
            dbg!(&challenge_set[i].operator_norm());
        }

        for i in 0..multiplicity {
            for j in 0..multiplicity {
                if i != j {
                    assert_ne!(challenge_set[i], challenge_set[j]);
                }
            }
        }
    }

    // TODO: Testing randomness of the distribution is needed.
}
