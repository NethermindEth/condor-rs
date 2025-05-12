use super::shake_sponge::Sponge;
use crate::ring::{rq::Rq, rq_vector::RqVector, zq::Zq};

pub struct LabradorTranscript<S: Sponge> {
    sponge: S,
    security_parameter: usize,
    rank: usize,
    multiplicity: usize,
    pub u1: RqVector,
    pub vector_p: Vec<Zq>,
    pub b_ct_aggr: RqVector,
    u2: RqVector,
}

impl<S: Sponge> LabradorTranscript<S> {
    pub fn new(sponge: S, security_parameter: usize, rank: usize, multiplicity: usize) -> Self {
        Self {
            sponge,
            security_parameter,
            rank,
            multiplicity,
            u1: RqVector::new(Vec::new()),
            vector_p: Vec::new(),
            b_ct_aggr: RqVector::new(Vec::new()),
            u2: RqVector::new(Vec::new()),
        }
    }

    pub fn absorb_u1(&mut self, u1: RqVector) {
        self.sponge.absorb_rq(u1.get_elements());
        self.u1 = u1;
    }

    pub fn absorb_vector_p(&mut self, p: Vec<Zq>) {
        self.sponge.absorb_zq(&p);
        self.vector_p = p;
    }

    pub fn absorb_vector_b_ct_aggr(&mut self, input: RqVector) {
        self.sponge.absorb_rq(input.get_elements());
        self.b_ct_aggr = input;
    }

    pub fn generate_vector_of_projection_matrices(&mut self) -> Vec<Vec<Vec<Zq>>> {
        // r vectors, each of length 256 * nD
        let row_size = 2 * self.security_parameter;
        let col_size = self.rank * Rq::DEGREE;

        (0..self.multiplicity)
            .map(|_| {
                let linear_projection_randomness = self.sponge.squeeze_zq(row_size * col_size);
                linear_projection_randomness
                    .chunks_exact(col_size)
                    .map(|chunk| {
                        chunk
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
                            .collect()
                    })
                    .collect()
            })
            .collect()
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

    pub fn generate_vector_omega(&mut self, number_of_vectors: usize) -> Vec<Vec<Zq>> {
        let elements = self
            .sponge
            .squeeze_zq(number_of_vectors * 2 * self.security_parameter);
        elements
            .chunks_exact(2 * self.security_parameter)
            .map(|chunk| chunk.into())
            .collect()
    }

    pub fn generate_rq_vector(&mut self, vector_length: usize) -> RqVector {
        RqVector::new(self.sponge.squeeze_rq(vector_length))
    }

    pub fn generate_ci(&mut self) {}
}

#[cfg(test)]
mod tests_generate_pi {
    use super::*;
    use crate::transcript::shake_sponge::ShakeSponge;

    #[test]
    fn test_projection_matrix_has_correct_size() {
        let (security_parameter, rank, multiplicity) = (128, 20, 9);
        let mut transcript = LabradorTranscript::new(
            ShakeSponge::default(),
            security_parameter,
            rank,
            multiplicity,
        );
        let result = transcript.generate_vector_of_projection_matrices();
        assert_eq!(result.len(), multiplicity); // number_of_project_matrices
        assert_eq!(result[0].len(), 2 * security_parameter);
        assert_eq!(result[0][0].len(), rank * Rq::DEGREE);
    }

    // Test the distribution of values in the random matrix
    #[test]
    #[allow(clippy::as_conversions)]
    fn test_projection_matrix_is_random() {
        let (security_parameter, rank, multiplicity) = (128, 1000, 1);
        let mut transcript = LabradorTranscript::new(
            ShakeSponge::default(),
            security_parameter,
            rank,
            multiplicity,
        );
        let projection_matrix_vector = transcript.generate_vector_of_projection_matrices();

        for projection_matrix in projection_matrix_vector {
            let mut counts = [0.0, 0.0, 0.0]; // -1, 0, 1
            for row in projection_matrix {
                for cell in row {
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
    use crate::transcript::shake_sponge::ShakeSponge;

    #[test]
    fn test_generate_vector_psi_has_correct_size() {
        let (security_parameter, rank, multiplicity, k_range, l_range) = (128, 20, 9, 20, 12);
        let mut transcript = LabradorTranscript::new(
            ShakeSponge::default(),
            security_parameter,
            rank,
            multiplicity,
        );
        let result = transcript.generate_vector_psi(k_range, l_range);
        assert_eq!(result.len(), k_range); // number_of_project_matrices
        assert_eq!(result[0].len(), l_range);
    }

    // TODO: Testing randomness of the distribution is needed.
}

#[cfg(test)]
mod test_generate_omega {
    use super::*;
    use crate::transcript::shake_sponge::ShakeSponge;

    #[test]
    fn test_generate_vector_omega_has_correct_size() {
        let (security_parameter, rank, multiplicity, k_range) = (128, 20, 9, 20);
        let mut transcript = LabradorTranscript::new(
            ShakeSponge::default(),
            security_parameter,
            rank,
            multiplicity,
        );
        let result = transcript.generate_vector_omega(k_range);
        assert_eq!(result.len(), k_range); // number_of_project_matrices
        assert_eq!(result[0].len(), 256);
    }

    // TODO: Testing randomness of the distribution is needed.
}

#[cfg(test)]
mod test_generate_rq {
    use super::*;
    use crate::transcript::shake_sponge::ShakeSponge;

    #[test]
    fn test_generate_rq_has_correct_size() {
        let (security_parameter, rank, multiplicity, k_range) = (128, 20, 9, 20);
        let mut transcript = LabradorTranscript::new(
            ShakeSponge::default(),
            security_parameter,
            rank,
            multiplicity,
        );
        let result = transcript.generate_rq_vector(k_range);
        assert_eq!(result.get_elements().len(), k_range); // number_of_project_matrices
    }

    // TODO: Testing randomness of the distribution is needed.
}
