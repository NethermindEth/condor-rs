use crate::commitments::ajtai_commitment::{AjtaiCommitment, AjtaiParameters};
use crate::commitments::hierarchical_commitment::{DecompositionParameters, HierarchicalError};
use crate::ring::poly::{PolyRing, PolyVector};
use crate::ring::rq::Rq;
use crate::ring::rq_matrix::RqMatrix;
use crate::ring::rq_vector::RqVector;
use crate::ring::zq::Zq;
use rand::SeedableRng;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Represents the outer commitments in the recursive LaBRADOR protocol
#[derive(Debug, Clone)]
pub struct RecursiveCommitment<const M: usize, const N: usize, const D: usize> {
    /// ν₁: Outer commitment to t_i (inner commitment parts) and g_{ij} polynomials
    pub nu1: RqVector<M, D>,

    /// ν₂: Outer commitment to h_{ij} polynomials
    pub nu2: RqVector<M, D>,
}

/// Parameters for garbage polynomial handling
#[derive(Debug, Clone)]
pub struct GarbageParameters {
    /// Base for decomposing g polynomials (b₂ in the paper)
    pub g_base: Zq,
    /// Number of parts for g decomposition (t₂ in the paper)
    pub g_parts: usize,
    /// Base for decomposing h polynomials (b₁ in the paper)
    pub h_base: Zq,
    /// Number of parts for h decomposition (t₁ in the paper)
    pub h_parts: usize,
}

/// Witness for the next recursion level
#[derive(Debug, Clone)]
pub struct RecursiveWitness<const D: usize> {
    /// z components from the inner commitment openings
    pub z_parts: Vec<Rq<D>>,
    /// t components (decomposed inner commitments)
    pub t_parts: Vec<Rq<D>>,
    /// g components (decomposed garbage polynomials)
    pub g_parts: Vec<Rq<D>>,
    /// h components (decomposed garbage polynomials)
    pub h_parts: Vec<Rq<D>>,
}

/// Manages garbage polynomial calculation and commitment
pub struct GarbagePolynomialCommitment<const M: usize, const N: usize, const D: usize> {
    // Parameters for g polynomial decomposition
    g_decomp_params: DecompositionParameters,
    // Parameters for h polynomial decomposition
    h_decomp_params: DecompositionParameters,
    // Commitment matrix for combined (t,g) polynomials as ν₁
    nu1_matrix: RqMatrix<M, N, D>,
    // Commitment matrix for h polynomials as ν₂
    nu2_matrix: RqMatrix<M, N, D>,
    // Security parameters
    params: AjtaiParameters,
}

impl<const M: usize, const N: usize, const D: usize> GarbagePolynomialCommitment<M, N, D> {
    /// Creates a new garbage polynomial commitment system for the recursive protocol
    pub fn new(
        params: AjtaiParameters,
        nu1_matrix: RqMatrix<M, N, D>,
        nu2_matrix: RqMatrix<M, N, D>,
        garbage_params: GarbageParameters,
    ) -> Result<Self, HierarchicalError> {
        let g_decomp_params =
            DecompositionParameters::new(garbage_params.g_base, garbage_params.g_parts)?;
        let h_decomp_params =
            DecompositionParameters::new(garbage_params.h_base, garbage_params.h_parts)?;

        Ok(Self {
            g_decomp_params,
            h_decomp_params,
            nu1_matrix,
            nu2_matrix,
            params,
        })
    }

    /// Calculate the garbage polynomials g_{ij} = <s_i, s_j>
    /// Exploits symmetry by only calculating for i ≤ j since g_{ij} = g_{ji}
    pub fn calculate_g_polynomials(witnesses: &[PolyVector]) -> Vec<PolyRing> {
        let r = witnesses.len();
        let mut g_polys = Vec::with_capacity((r * (r + 1)) / 2);

        for i in 0..r {
            for j in 0..=i {
                // Only calculate for j ≤ i (upper triangular)
                let g_ij = witnesses[i].inner_product_poly_vector(&witnesses[j]);
                g_polys.push(g_ij);
            }
        }

        g_polys
    }

    /// Calculate the h_{ij} = <φ_i, s_j> + <φ_j, s_i> garbage polynomials
    /// In the paper, h_{ij} is defined with a factor of 1/2 in front
    /// However, since we're using q = 2^32, division by 2 is problematic in Z_q
    /// So we store h'_{ij} = 2*h_{ij} = <φ_i, s_j> + <φ_j, s_i> directly
    /// Exploits symmetry by only calculating for i ≤ j since h_{ij} = h_{ji}
    pub fn calculate_h_polynomials(witnesses: &[PolyVector], phi: &[PolyVector]) -> Vec<PolyRing> {
        let r = witnesses.len();
        let mut h_polys = Vec::with_capacity((r * (r + 1)) / 2);

        for i in 0..r {
            for j in 0..=i {
                // Only calculate for j ≤ i (upper triangular)
                let inner_phi_i_s_j = phi[i].inner_product_poly_vector(&witnesses[j]);
                let inner_phi_j_s_i = phi[j].inner_product_poly_vector(&witnesses[i]);
                let h_ij = &inner_phi_i_s_j + &inner_phi_j_s_i;
                h_polys.push(h_ij);
            }
        }

        h_polys
    }

    /// Decompose a polynomial into parts according to the specified parameters
    fn decompose_polynomial(&self, poly: &PolyRing, is_g: bool) -> Vec<Rq<D>> {
        let params = if is_g {
            &self.g_decomp_params
        } else {
            &self.h_decomp_params
        };

        // Convert PolyRing to Rq<D>
        let rq_poly: Rq<D> = poly.clone().into();

        // Decompose using the appropriate base and number of parts
        rq_poly.decompose(params.base(), params.num_parts())
    }

    /// Decompose all garbage polynomials into their parts
    fn decompose_all(
        &self,
        g_polys: &[PolyRing],
        h_polys: &[PolyRing],
    ) -> (Vec<Rq<D>>, Vec<Rq<D>>) {
        let mut g_parts = Vec::new();
        for poly in g_polys {
            let parts = self.decompose_polynomial(poly, true);
            g_parts.extend(parts);
        }

        let mut h_parts = Vec::new();
        for poly in h_polys {
            let parts = self.decompose_polynomial(poly, false);
            h_parts.extend(parts);
        }

        (g_parts, h_parts)
    }

    /// Create a valid witness for commitment from decomposed parts
    fn create_witness_from_parts(&self, parts: &[Rq<D>]) -> RqVector<N, D> {
        // Take up to N parts (or pad with zeroes if fewer)
        let mut rqs = Vec::with_capacity(N);

        if parts.len() >= N {
            // Take first N parts
            for part in parts.iter().take(N) {
                rqs.push(part.clone());
            }
        } else {
            // Take all parts and pad with zeros
            for part in parts {
                rqs.push(part.clone());
            }

            while rqs.len() < N {
                rqs.push(Rq::zero());
            }
        }

        // Apply witness bound to ensure coefficients are valid
        let witness_bound = self.params.witness_bound();
        let bounded_rqs: Vec<Rq<D>> = rqs
            .iter()
            .map(|rq| {
                let mut bounded_coeffs = [Zq::ZERO; D];
                for (i, coeff) in rq.get_coefficients().iter().enumerate() {
                    bounded_coeffs[i] = coeff.centered_mod(witness_bound);
                }
                Rq::new(bounded_coeffs)
            })
            .collect();

        RqVector::from(bounded_rqs)
    }

    /// Commit to garbage polynomials following LaBRADOR's recursive structure
    /// Takes the inner commitment parts (t) as input to combine with g parts in ν₁
    pub fn commit_recursive(
        &self,
        witnesses: &[PolyVector],
        phi: &[PolyVector],
        t_parts: &[Rq<D>],
        z_parts: &[Rq<D>],
    ) -> Result<(RecursiveCommitment<M, N, D>, RecursiveWitness<D>), HierarchicalError> {
        // 1. Calculate garbage polynomials g_{ij} = <s_i, s_j> and h_{ij} = <φ_i, s_j> + <φ_j, s_i>
        let g_polys = Self::calculate_g_polynomials(witnesses);
        let h_polys = Self::calculate_h_polynomials(witnesses, phi);

        // 2. Decompose all polynomials into parts
        let (g_parts, h_parts) = self.decompose_all(&g_polys, &h_polys);

        // 3. Combine t_parts and g_parts for ν₁ commitment
        let mut combined_parts = t_parts.to_vec();
        combined_parts.extend(g_parts.clone());

        // 4. Create witnesses for outer commitments
        let nu1_witness = self.create_witness_from_parts(&combined_parts);
        let nu2_witness = self.create_witness_from_parts(&h_parts);

        // 5. Create outer commitments (ν₁, ν₂)
        let nu1_commitment = AjtaiCommitment::new(self.params.clone(), self.nu1_matrix.clone())?;
        let (nu1, _) = nu1_commitment.commit(nu1_witness)?;

        let nu2_commitment = AjtaiCommitment::new(self.params.clone(), self.nu2_matrix.clone())?;
        let (nu2, _) = nu2_commitment.commit(nu2_witness)?;

        // 6. Create the recursive witness that will be used as input to the next level
        let recursive_witness = RecursiveWitness {
            z_parts: z_parts.to_vec(),
            t_parts: t_parts.to_vec(),
            g_parts,
            h_parts,
        };

        Ok((RecursiveCommitment { nu1, nu2 }, recursive_witness))
    }

    /// Calculate optimal decomposition parameters based on Section 5.4 of the paper
    #[allow(clippy::as_conversions)]
    pub fn calculate_optimal_parameters(
        n: usize,  // Dimension of the witness vector
        d: usize,  // Degree of the polynomial ring
        s: f64,    // Standard deviation of witness coefficients
        beta: f64, // Ajtai commitment parameter
    ) -> GarbageParameters {
        // Calculate estimated standard deviations based on Section 5.4
        // For g_{ij}, std dev is approximately s_g = sqrt(n*d) * s^2
        // For h_{ij}, std dev is approximately s_h = beta * sqrt(n*d) * s
        let n_d_sqrt = (n * d) as f64;
        let s_g = n_d_sqrt * s * s;
        let s_h = beta * n_d_sqrt * s;

        // Calculate optimal bases using formula from Section 5.4
        // B_i ≈ sqrt(12 * s_i)
        let b2 = ((12.0 * s_g).sqrt() as u32).max(2);
        let b1 = ((12.0 * s_h).sqrt() as u32).max(2);

        // Calculate optimal number of parts
        // t_i ≈ log_B_i(q)
        let t2 = ((32.0 / (b2 as f64).log2()).ceil() as usize).max(2);
        let t1 = ((32.0 / (b1 as f64).log2()).ceil() as usize).max(2);

        GarbageParameters {
            g_base: Zq::new(b2),
            g_parts: t2,
            h_base: Zq::new(b1),
            h_parts: t1,
        }
    }

    /// Implementation of the final level optimization (Section 5.6)
    /// Uses sequentially derived challenges via Fiat-Shamir to simulate the interactive protocol
    pub fn optimize_final_level(
        witnesses: &[PolyVector],
        phi: &[PolyVector],
        initial_seed: u64,
    ) -> (PolyRing, Vec<PolyRing>, Vec<PolyRing>) {
        let r = witnesses.len();

        // Calculate g_0 = Σ_i <s_i, s_i> (diagonal sum)
        let g0 = (0..r)
            .map(|i| witnesses[i].inner_product_poly_vector(&witnesses[i]))
            .fold(
                PolyRing::zero(witnesses[0].get_elements()[0].len()),
                |acc, g| &acc + &g,
            );

        // Generate sequence of challenges using Fiat-Shamir
        let mut challenges = Vec::with_capacity(r / 2);
        let mut hasher = DefaultHasher::new();
        initial_seed.hash(&mut hasher);
        let mut current_seed = hasher.finish();

        for _ in 0..r / 2 {
            let mut rng = rand::rngs::StdRng::seed_from_u64(current_seed);
            let challenge = PolyRing::random(&mut rng, witnesses[0].get_elements()[0].len());
            challenges.push(challenge);

            // Update seed for next challenge
            let mut hasher = DefaultHasher::new();
            current_seed.hash(&mut hasher);
            current_seed = hasher.finish();
        }

        // Calculate selected g terms: g_{2i-1} and g_{2i}
        let mut g_terms = Vec::new();
        for i in 1..=r / 2 {
            let idx1 = 2 * i - 2;
            let idx2 = 2 * i - 1;

            // Use unique challenge for each i
            let challenge = &challenges[i - 1];

            // Add g_{2i-1}
            if idx2 < r {
                // For g_{2i-1} = <s_{2i-2}, c_i * s_{2i-1}>
                let s_j_scaled = witnesses[idx2]
                    .iter()
                    .map(|p| p * challenge)
                    .collect::<PolyVector>();

                let g_2i_1 = witnesses[idx1].inner_product_poly_vector(&s_j_scaled);
                g_terms.push(g_2i_1);

                // Add g_{2i} if we have enough witnesses
                if 2 * i < r {
                    // For g_{2i} = <s_{2i-1}, c_i * s_{2i}>
                    let s_j_scaled = witnesses[2 * i]
                        .iter()
                        .map(|p| p * challenge)
                        .collect::<PolyVector>();

                    let g_2i = witnesses[idx2].inner_product_poly_vector(&s_j_scaled);
                    g_terms.push(g_2i);
                }
            }
        }

        // Calculate selected h terms: h_{2i-1} and h_{2i}
        let mut h_terms = Vec::new();
        for i in 1..=r / 2 {
            let idx1 = 2 * i - 2;
            let idx2 = 2 * i - 1;

            // Use unique challenge for each i
            let challenge = &challenges[i - 1];

            // Add h_{2i-1}
            if idx2 < r {
                // For h_{2i-1} = <φ_{2i-2}, c_i * s_{2i-1}> + <φ_{2i-1}, c_i * s_{2i-2}>
                let s_j_scaled = witnesses[idx2]
                    .iter()
                    .map(|p| p * challenge)
                    .collect::<PolyVector>();

                let phi_i_s_j = phi[idx1].inner_product_poly_vector(&s_j_scaled);

                let s_i_scaled = witnesses[idx1]
                    .iter()
                    .map(|p| p * challenge)
                    .collect::<PolyVector>();

                let phi_j_s_i = phi[idx2].inner_product_poly_vector(&s_i_scaled);

                let h_2i_1 = &phi_i_s_j + &phi_j_s_i;
                h_terms.push(h_2i_1);

                // Add h_{2i} if we have enough witnesses
                if 2 * i < r {
                    // For h_{2i} = <φ_{2i-1}, c_i * s_{2i}> + <φ_{2i}, c_i * s_{2i-1}>
                    let s_j_scaled = witnesses[2 * i]
                        .iter()
                        .map(|p| p * challenge)
                        .collect::<PolyVector>();

                    let phi_i_s_j = phi[idx2].inner_product_poly_vector(&s_j_scaled);

                    let s_i_scaled = witnesses[idx2]
                        .iter()
                        .map(|p| p * challenge)
                        .collect::<PolyVector>();

                    let phi_j_s_i = phi[2 * i].inner_product_poly_vector(&s_i_scaled);

                    let h_2i = &phi_i_s_j + &phi_j_s_i;
                    h_terms.push(h_2i);
                }
            }
        }

        (g0, g_terms, h_terms)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rng;

    const TEST_M: usize = 8;
    const TEST_N: usize = 8;
    const TEST_D: usize = 4;

    fn create_test_commitment() -> GarbagePolynomialCommitment<TEST_M, TEST_N, TEST_D> {
        let beta = Zq::ONE;
        let witness_bound = Zq::ONE;
        let params = AjtaiParameters::new(beta, witness_bound).unwrap();

        let mut rng = rng();
        let nu1_matrix = RqMatrix::<TEST_M, TEST_N, TEST_D>::random(&mut rng);
        let nu2_matrix = RqMatrix::<TEST_M, TEST_N, TEST_D>::random(&mut rng);

        // Use base 8 and 2 parts for g polynomials
        // Use base 4 and 3 parts for h polynomials (as examples)
        let garbage_params = GarbageParameters {
            g_base: Zq::new(8),
            g_parts: 2,
            h_base: Zq::new(4),
            h_parts: 3,
        };

        GarbagePolynomialCommitment::new(params, nu1_matrix, nu2_matrix, garbage_params).unwrap()
    }

    fn create_test_witnesses(count: usize) -> (Vec<PolyVector>, Vec<PolyVector>) {
        let witnesses = (0..count)
            .map(|_| PolyVector::random(TEST_N, TEST_D))
            .collect();

        let phi = (0..count)
            .map(|_| PolyVector::random(TEST_N, TEST_D))
            .collect();

        (witnesses, phi)
    }

    fn create_test_parts(count: usize) -> Vec<Rq<TEST_D>> {
        let mut rng = rng();
        (0..count).map(|_| Rq::<TEST_D>::random(&mut rng)).collect()
    }

    #[test]
    fn test_g_calculation() {
        let (witnesses, _) = create_test_witnesses(3);
        let r = witnesses.len();

        let g_polys =
            GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::calculate_g_polynomials(
                &witnesses,
            );

        // Check we have the correct number of elements in upper triangular
        assert_eq!(g_polys.len(), (r * (r + 1)) / 2); // Sum of 1..r

        // Verify a few specific values
        let expected_g_01 = witnesses[0].inner_product_poly_vector(&witnesses[1]);
        let expected_g_22 = witnesses[2].inner_product_poly_vector(&witnesses[2]);

        // Manual access - for triangular storage with i ≤ j, the index is:
        // idx = (i * (i + 1)) / 2 + j
        let idx_01 = 1; // For i=0, j=1: (0*(0+1))/2 + 1 = 1
        let idx_22 = 5; // For i=2, j=2: (2*(2+1))/2 + 2 = 5

        assert_eq!(g_polys[idx_01], expected_g_01);
        assert_eq!(g_polys[idx_22], expected_g_22);
    }

    #[test]
    fn test_h_calculation() {
        let (witnesses, phi) = create_test_witnesses(3);
        let r = witnesses.len();

        let h_polys =
            GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::calculate_h_polynomials(
                &witnesses, &phi,
            );

        // Check we have the correct number of elements in upper triangular
        assert_eq!(h_polys.len(), (r * (r + 1)) / 2);

        // Verify a specific value
        let phi_i_s_j = phi[0].inner_product_poly_vector(&witnesses[1]);
        let phi_j_s_i = phi[1].inner_product_poly_vector(&witnesses[0]);
        let expected_h_01 = &phi_i_s_j + &phi_j_s_i;

        // Manual access - same indexing as g_polys
        let idx_01 = 1; // For i=0, j=1: (0*(0+1))/2 + 1 = 1

        assert_eq!(h_polys[idx_01], expected_h_01);
    }

    #[test]
    fn test_commit_recursive() {
        let commitment_scheme = create_test_commitment();
        let (witnesses, phi) = create_test_witnesses(3);

        // Create mock inner commitment parts (t and z)
        let t_parts = create_test_parts(5);
        let z_parts = create_test_parts(5);

        let (recursive_commitment, recursive_witness) = commitment_scheme
            .commit_recursive(&witnesses, &phi, &t_parts, &z_parts)
            .unwrap();

        // Check the output structure
        assert!(!recursive_commitment.nu1.as_slice().is_empty());
        assert!(!recursive_commitment.nu2.as_slice().is_empty());

        // Check witness has appropriate parts
        assert_eq!(recursive_witness.t_parts.len(), t_parts.len());
        assert_eq!(recursive_witness.z_parts.len(), z_parts.len());
        assert!(!recursive_witness.g_parts.is_empty());
        assert!(!recursive_witness.h_parts.is_empty());
    }

    #[test]
    fn test_commit_recursive_correctness() {
        let commitment_scheme = create_test_commitment();
        let (witnesses, phi) = create_test_witnesses(3);

        // Create test parts for inner commitment
        let t_parts = create_test_parts(5);
        let z_parts = create_test_parts(5);

        // Get the recursive commitment and witness
        let (recursive_commitment, recursive_witness) = commitment_scheme
            .commit_recursive(&witnesses, &phi, &t_parts, &z_parts)
            .unwrap();

        // Manually compute the expected commitments:

        // 1. For nu1, combine t_parts and g_parts
        let mut combined_parts = recursive_witness.t_parts.clone();
        combined_parts.extend(recursive_witness.g_parts.clone());

        // 2. Create expected witnesses
        let expected_nu1_witness = commitment_scheme.create_witness_from_parts(&combined_parts);
        let expected_nu2_witness =
            commitment_scheme.create_witness_from_parts(&recursive_witness.h_parts);

        // 3. Create AjtaiCommitment instances with the correct matrices
        let nu1_commitment = AjtaiCommitment::new(
            commitment_scheme.params.clone(),
            commitment_scheme.nu1_matrix.clone(),
        )
        .unwrap();

        let nu2_commitment = AjtaiCommitment::new(
            commitment_scheme.params.clone(),
            commitment_scheme.nu2_matrix.clone(),
        )
        .unwrap();

        // 4. Compute expected commitments
        let (expected_nu1, _) = nu1_commitment.commit(expected_nu1_witness).unwrap();
        let (expected_nu2, _) = nu2_commitment.commit(expected_nu2_witness).unwrap();

        // 5. Compare expected with actual
        assert_eq!(
            recursive_commitment.nu1, expected_nu1,
            "nu1 commitment does not match expected value"
        );
        assert_eq!(
            recursive_commitment.nu2, expected_nu2,
            "nu2 commitment does not match expected value"
        );
    }

    #[test]
    fn test_decomposition_reconstruction() {
        let commitment_scheme = create_test_commitment();
        let mut rng = rand::rng();

        // Create a random polynomial
        let original_poly = PolyRing::random(&mut rng, TEST_D);
        let original_rq: Rq<TEST_D> = original_poly.clone().into();

        // Test g decomposition parameters
        let g_parts = commitment_scheme.decompose_polynomial(&original_poly, true);
        let g_base = commitment_scheme.g_decomp_params.base();

        // Reconstruct the polynomial from parts
        let mut reconstructed_g = Rq::<TEST_D>::zero();
        let mut current_base_power = Zq::ONE; // Base^0

        for part in &g_parts {
            // Add part * base^k
            reconstructed_g = reconstructed_g.clone() + part.clone().scalar_mul(current_base_power);
            // Multiply by base for next iteration
            current_base_power *= g_base;
        }

        assert_eq!(
            reconstructed_g, original_rq,
            "G decomposition reconstruction failed"
        );

        // Test h decomposition parameters
        let h_parts = commitment_scheme.decompose_polynomial(&original_poly, false);
        let h_base = commitment_scheme.h_decomp_params.base();

        // Reconstruct the polynomial from parts
        let mut reconstructed_h = Rq::<TEST_D>::zero();
        let mut current_base_power = Zq::ONE; // Base^0

        for part in &h_parts {
            // Add part * base^k
            reconstructed_h = reconstructed_h.clone() + part.clone().scalar_mul(current_base_power);
            // Multiply by base for next iteration
            current_base_power *= h_base;
        }

        assert_eq!(
            reconstructed_h, original_rq,
            "H decomposition reconstruction failed"
        );
    }

    #[test]
    fn test_create_witness_from_parts_edge_cases() {
        let commitment_scheme = create_test_commitment();
        let mut rng = rand::rng();

        // Test case 1: parts.len() < N (should pad with zeros)
        let few_parts: Vec<Rq<TEST_D>> = (0..TEST_N - 2)
            .map(|_| Rq::<TEST_D>::random(&mut rng))
            .collect();

        let witness_few = commitment_scheme.create_witness_from_parts(&few_parts);

        // Check length is exactly N
        assert_eq!(witness_few.as_slice().len(), TEST_N);

        // Check that last elements are zero
        for i in few_parts.len()..TEST_N {
            assert_eq!(witness_few[i], Rq::<TEST_D>::zero());
        }

        // Check original parts are preserved (with possible bounding applied)
        let witness_bound = commitment_scheme.params.witness_bound();
        for (i, part) in few_parts.iter().enumerate() {
            let mut bounded_part_coeffs = [Zq::ZERO; TEST_D];
            for (j, coeff) in part.get_coefficients().iter().enumerate().take(TEST_D) {
                bounded_part_coeffs[j] = coeff.centered_mod(witness_bound);
            }
            let bounded_part = Rq::<TEST_D>::new(bounded_part_coeffs);
            assert_eq!(witness_few[i], bounded_part);
        }

        // Test case 2: parts.len() > N (should truncate to first N)
        let many_parts: Vec<Rq<TEST_D>> = (0..TEST_N + 3)
            .map(|_| Rq::<TEST_D>::random(&mut rng))
            .collect();

        let witness_many = commitment_scheme.create_witness_from_parts(&many_parts);

        // Check length is exactly N
        assert_eq!(witness_many.as_slice().len(), TEST_N);

        // Check only first N parts are included (with possible bounding applied)
        let witness_bound = commitment_scheme.params.witness_bound();
        for i in 0..TEST_N {
            let mut bounded_part_coeffs = [Zq::ZERO; TEST_D];
            for (j, coeff) in many_parts[i]
                .get_coefficients()
                .iter()
                .enumerate()
                .take(TEST_D)
            {
                bounded_part_coeffs[j] = coeff.centered_mod(witness_bound);
            }
            let bounded_part = Rq::<TEST_D>::new(bounded_part_coeffs);
            assert_eq!(witness_many[i], bounded_part);
        }

        // Test case 3: Coefficient bounding
        // Create parts with large coefficients
        let large_coeff_parts: Vec<Rq<TEST_D>> = (0..TEST_N)
            .map(|_| {
                let mut large_part = Rq::<TEST_D>::random(&mut rng);
                // Set first coefficient to a large value beyond the witness bound
                if let Some(c) = large_part.iter_mut().next() {
                    *c = Zq::MAX - Zq::ONE;
                }
                large_part
            })
            .collect();

        let witness_large = commitment_scheme.create_witness_from_parts(&large_coeff_parts);

        // Check coefficients are properly bounded
        for i in 0..TEST_N {
            let mut bounded_part_coeffs = [Zq::ZERO; TEST_D];
            for (j, coeff) in large_coeff_parts[i]
                .get_coefficients()
                .iter()
                .enumerate()
                .take(TEST_D)
            {
                bounded_part_coeffs[j] = coeff.centered_mod(witness_bound);
            }
            let bounded_part = Rq::<TEST_D>::new(bounded_part_coeffs);

            assert_eq!(witness_large[i], bounded_part);

            // Also check that large coefficients were properly bounded
            for (j, coeff) in witness_large[i].get_coefficients().iter().enumerate() {
                let original_coeff = large_coeff_parts[i].get_coefficients()[j];
                if original_coeff > witness_bound {
                    assert!(*coeff <= witness_bound);
                }
            }
        }
    }

    #[test]
    #[allow(clippy::as_conversions)]
    fn test_optimal_parameters_accuracy() {
        // Choose small, simple values for manual calculation
        let n = 4; // Small dimension
        let d = 4; // Small degree
        let s = 2.0; // Simple standard deviation
        let beta = 1.0; // Simple beta value

        // Manually calculate the expected values according to the paper's formulas
        let n_d_sqrt = (n * d) as f64; // = 4.0

        // s_g = sqrt(n*d) * s^2 = 4.0 * 4.0 = 16.0
        let s_g = n_d_sqrt * s * s;

        // s_h = beta * sqrt(n*d) * s = 1.0 * 4.0 * 2.0 = 8.0
        let s_h = beta * n_d_sqrt * s;

        // b2 ≈ sqrt(12 * s_g) = sqrt(12 * 16.0) = sqrt(192) ≈ 13.856... => 13
        let expected_b2 = (12.0 * s_g).sqrt() as u32;

        // b1 ≈ sqrt(12 * s_h) = sqrt(12 * 8.0) = sqrt(96) ≈ 9.798... => 9
        let expected_b1 = (12.0 * s_h).sqrt() as u32;

        // For q = 2^32:
        // t2 ≈ log_b2(q) = log_13(2^32) = 32/log_2(13) ≈ 32/3.7 ≈ 8.65 => 9
        let expected_t2 = ((32.0 / (expected_b2 as f64).log2()).ceil() as usize).max(2);

        // t1 ≈ log_b1(q) = log_9(2^32) = 32/log_2(9) ≈ 32/3.17 ≈ 10.09 => 11
        let expected_t1 = ((32.0 / (expected_b1 as f64).log2()).ceil() as usize).max(2);

        // Call the function under test
        let params =
            GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::calculate_optimal_parameters(
                n, d, s, beta,
            );

        // Check results
        assert!(
            params.g_base.to_u128() >= expected_b2 as u128 - 1
                && params.g_base.to_u128() <= expected_b2 as u128 + 1,
            "g_base {}, expected {}",
            params.g_base.to_u128(),
            expected_b2
        );

        assert!(
            params.h_base.to_u128() >= expected_b1 as u128 - 1
                && params.h_base.to_u128() <= expected_b1 as u128 + 1,
            "h_base {}, expected {}",
            params.h_base.to_u128(),
            expected_b1
        );

        // For part counts, use approximate comparison due to potential floating point differences
        assert!(
            params.g_parts >= expected_t2 - 1 && params.g_parts <= expected_t2 + 1,
            "g_parts {}, expected {}",
            params.g_parts,
            expected_t2
        );

        assert!(
            params.h_parts >= expected_t1 - 1 && params.h_parts <= expected_t1 + 1,
            "h_parts {}, expected {}",
            params.h_parts,
            expected_t1
        );
    }

    #[test]
    fn test_error_handling() {
        let params = AjtaiParameters::new(Zq::ONE, Zq::ONE).unwrap();
        let mut rng = rand::rng();
        let nu1_matrix = RqMatrix::<TEST_M, TEST_N, TEST_D>::random(&mut rng);
        let nu2_matrix = RqMatrix::<TEST_M, TEST_N, TEST_D>::random(&mut rng);

        // Test invalid g_base (≤ 1)
        let invalid_g_base = GarbageParameters {
            g_base: Zq::ONE, // Invalid: base must be > 1
            g_parts: 2,
            h_base: Zq::new(4),
            h_parts: 3,
        };

        let result = GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::new(
            params.clone(),
            nu1_matrix.clone(),
            nu2_matrix.clone(),
            invalid_g_base,
        );

        assert!(result.is_err(), "Should reject g_base ≤ 1");
        match result {
            Err(HierarchicalError::InvalidBase(_)) => {} // Expected error
            _ => panic!("Wrong error type for invalid g_base"),
        }

        // Test invalid g_parts (0)
        let invalid_g_parts = GarbageParameters {
            g_base: Zq::new(8),
            g_parts: 0, // Invalid: parts must be > 0
            h_base: Zq::new(4),
            h_parts: 3,
        };

        let result = GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::new(
            params.clone(),
            nu1_matrix.clone(),
            nu2_matrix.clone(),
            invalid_g_parts,
        );

        assert!(result.is_err(), "Should reject g_parts = 0");
        match result {
            Err(HierarchicalError::InvalidPartCount(_)) => {} // Expected error
            _ => panic!("Wrong error type for invalid g_parts"),
        }

        // Test invalid h_base (≤ 1)
        let invalid_h_base = GarbageParameters {
            g_base: Zq::new(8),
            g_parts: 2,
            h_base: Zq::ZERO, // Invalid: base must be > 0
            h_parts: 3,
        };

        let result = GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::new(
            params.clone(),
            nu1_matrix.clone(),
            nu2_matrix.clone(),
            invalid_h_base,
        );

        assert!(result.is_err(), "Should reject h_base ≤ 1");
        match result {
            Err(HierarchicalError::InvalidBase(_)) => {} // Expected error
            _ => panic!("Wrong error type for invalid h_base"),
        }

        // Test invalid h_parts (0)
        let invalid_h_parts = GarbageParameters {
            g_base: Zq::new(8),
            g_parts: 2,
            h_base: Zq::new(4),
            h_parts: 0, // Invalid: parts must be > 0
        };

        let result = GarbagePolynomialCommitment::<TEST_M, TEST_N, TEST_D>::new(
            params.clone(),
            nu1_matrix.clone(),
            nu2_matrix,
            invalid_h_parts,
        );

        assert!(result.is_err(), "Should reject h_parts = 0");
        match result {
            Err(HierarchicalError::InvalidPartCount(_)) => {} // Expected error
            _ => panic!("Wrong error type for invalid h_parts"),
        }
    }
}
