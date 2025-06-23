//! Polynomial ring **R = Z_q[X]/(X^d + 1)** with `d = 64` and `q = u32::MAX`.

use crate::ring::zq::Zq;
use crate::ring::Norms;
use core::ops::{Add, Mul, Sub};
use rand::distr::{Distribution, Uniform};
use rand::{CryptoRng, Rng};
use rustfft::num_complex::Complex;
use rustfft::FftPlanner;

/// This module provides implementations for various operations
/// in the polynomial ring R = Z_q\[X\] / (X^d + 1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rq {
    coeffs: [Zq; Self::DEGREE],
}

impl Rq {
    pub const DEGREE: usize = 64;
    /// Constructor for the polynomial ring
    pub const fn new(coeffs: [Zq; Self::DEGREE]) -> Self {
        Self { coeffs }
    }

    /// Generate zero polynomial
    pub const fn zero() -> Self {
        Self {
            coeffs: [Zq::ZERO; Self::DEGREE],
        }
    }

    pub fn into_coeffs(self) -> [Zq; Self::DEGREE] {
        self.coeffs
    }

    /// Get the coefficients as a vector
    pub fn coeffs(&self) -> &[Zq; Self::DEGREE] {
        &self.coeffs
    }

    /// Random coefficients in `(0, Zq::MAX)`.
    pub fn random<R: Rng + CryptoRng>(rng: &mut R) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::NEG_ONE).unwrap();
        let mut coeffs = [Zq::ZERO; Self::DEGREE];
        coeffs.iter_mut().for_each(|c| *c = uniform.sample(rng));
        Self { coeffs }
    }

    /// Random coefficients in `(-bound, bound)`.
    pub fn random_with_bound<R: Rng + CryptoRng>(rng: &mut R, bound: u32) -> Self {
        let uniform = Uniform::new_inclusive(Zq::ZERO, Zq::new(bound)).unwrap();
        let mut coeffs = [Zq::ZERO; Self::DEGREE];
        coeffs.iter_mut().for_each(|c| {
            *c = if rng.random_bool(0.5) {
                -uniform.sample(rng)
            } else {
                uniform.sample(rng)
            }
        });
        Self { coeffs }
    }

    /// Decomposes a polynomial into base-B representation:
    /// p = p⁽⁰⁾ + p⁽¹⁾·B + p⁽²⁾·B² + ... + p⁽ᵗ⁻¹⁾·B^(t-1)
    /// Where each p⁽ⁱ⁾ has small coefficients, using centered representatives
    pub fn decompose(&self, base: Zq, num_parts: usize) -> Vec<Rq> {
        let mut result = vec![Rq::zero(); num_parts];
        self.coeffs.iter().enumerate().for_each(|(index, coeff)| {
            coeff
                .decompose(base, num_parts)
                .into_iter()
                .zip(result.iter_mut())
                .for_each(|(decomposed_coeff, decomposed_vec)| {
                    decomposed_vec.coeffs[index] = decomposed_coeff;
                });
        });
        result
    }

    /// Compute the conjugate automorphism \sigma_{-1} of vector based on B) Constraints..., Page 21.
    /// σ_{-1}(a_0, a_1, …, a_{n-1}) = (a_0, −a_{n-1}, −a_{n-2}, …, −a_1)
    pub fn conjugate_automorphism(&self) -> Self {
        let original_coefficients = self.coeffs();

        let mut conjugated_coeffs = [Zq::ZERO; Rq::DEGREE];
        conjugated_coeffs[0] = original_coefficients[0];

        for (conj_coeff, original_coeff) in conjugated_coeffs
            .iter_mut()
            .skip(1)
            .zip(original_coefficients.iter().rev())
        {
            *conj_coeff = original_coeff * Zq::NEG_ONE;
        }

        Self::new(conjugated_coeffs)
    }

    /// Compute the operator norm of a polynomial given its coefficients.
    /// The operator norm is defined as the maximum magnitude of the DFT (eigenvalues)
    /// of the coefficient vector.
    ///
    /// Note that: The operator norm only affects the coefficients of the random PolyRings generated from the challenge space.
    /// Prover and Verifier will not do the operator norm check, because random PolyRings are determined after generation.
    /// Both party will have access to the same PolyRings through transcript,
    #[allow(clippy::as_conversions)]
    pub fn operator_norm(&self) -> f64 {
        let coeffs = self.coeffs();
        let n = coeffs.len();
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);

        // Convert coefficients into complex numbers (with zero imaginary parts)
        let mut buffer: Vec<Complex<f64>> = coeffs
            .iter()
            .map(|&x| {
                let half = Zq::NEG_ONE.div_floor_by(2);
                let converted_value = if x > half {
                    x.to_u128() as f64 - Zq::NEG_ONE.to_u128() as f64 - 1.0
                } else {
                    x.to_u128() as f64
                };
                Complex {
                    re: converted_value,
                    im: 0.0,
                }
            })
            .collect();

        // Compute the FFT (this gives the eigenvalues of the circulant matrix)
        fft.process(&mut buffer);

        // Return the maximum absolute value (norm) among the eigenvalues
        buffer
            .iter()
            .map(|c| c.norm())
            .fold(0.0, |max, x| max.max(x))
    }
}

impl Add<&Rq> for &Rq {
    type Output = Rq;
    /// Add two polynomials
    fn add(self, other: &Rq) -> Rq {
        let mut coeffs = [Zq::ZERO; Rq::DEGREE];
        for (r, (a, b)) in coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a + *b;
        }
        Rq::new(coeffs)
    }
}

impl Sub<&Rq> for &Rq {
    type Output = Rq;
    /// Add two polynomials
    fn sub(self, other: &Rq) -> Rq {
        let mut coeffs = [Zq::ZERO; Rq::DEGREE];
        for (r, (a, b)) in coeffs
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = *a - *b;
        }
        Rq::new(coeffs)
    }
}

impl Mul<&Rq> for &Rq {
    type Output = Rq;
    /// Polynomial multiplication modulo x^D + 1
    fn mul(self, other: &Rq) -> Rq {
        let mut result = [Zq::ZERO; Rq::DEGREE];
        let mut out_of_field = [Zq::ZERO; Rq::DEGREE];
        for (i, &self_coeff) in self.coeffs.iter().enumerate() {
            for (j, &other_coeff) in other.coeffs.iter().enumerate() {
                if i + j < Rq::DEGREE {
                    result[i + j] += self_coeff * other_coeff;
                } else {
                    out_of_field[(i + j) % Rq::DEGREE] += self_coeff * other_coeff;
                }
            }
        }
        // Process excess terms with sign adjustment
        for i in (0..Rq::DEGREE).rev() {
            result[i] -= out_of_field[i];
        }
        Rq::new(result)
    }
}

impl Mul<&Zq> for &Rq {
    type Output = Rq;
    /// Scalar multiplication of a polynomial
    fn mul(self, other: &Zq) -> Rq {
        let mut copied_coeffs = self.coeffs;
        for elem in copied_coeffs.iter_mut() {
            *elem *= *other;
        }
        Rq::new(copied_coeffs)
    }
}

impl Norms for Rq {
    type NormType = u128;

    #[allow(clippy::as_conversions)]
    fn l2_norm_squared(&self) -> Self::NormType {
        self.coeffs.l2_norm_squared()
    }

    fn linf_norm(&self) -> Self::NormType {
        self.coeffs.linf_norm()
    }
}

#[cfg(test)]
mod norm_tests {
    use super::tests::generate_rq_from_zq_vector;
    use crate::ring::{rq::Rq, zq::Zq, Norms};

    // Test the square of the norm
    #[test]
    fn test_l2_norm() {
        let poly1 = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::NEG_ONE]);
        let poly2 = generate_rq_from_zq_vector(vec![Zq::ZERO, Zq::ZERO, Zq::new(5), Zq::ONE]);
        let poly3 = generate_rq_from_zq_vector(vec![Zq::new(5), Zq::ONE, -Zq::new(6), Zq::ZERO]);
        let poly4 = Rq::zero();

        assert_eq!(poly1.l2_norm_squared(), 27);
        assert_eq!(poly2.l2_norm_squared(), 26);
        assert_eq!(poly3.l2_norm_squared(), 62);
        assert_eq!(poly4.l2_norm_squared(), 0);
    }

    #[test]
    fn test_linf_norm() {
        let poly1 = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::NEG_ONE]);
        let poly2 = generate_rq_from_zq_vector(vec![Zq::ZERO, Zq::ZERO, -Zq::new(5), Zq::ONE]);
        let poly3 = generate_rq_from_zq_vector(vec![Zq::new(5), Zq::ONE, -Zq::new(6), Zq::ZERO]);
        let poly4 = Rq::zero();

        assert_eq!(poly1.linf_norm(), 5);
        assert_eq!(poly2.linf_norm(), 5);
        assert_eq!(poly3.linf_norm(), 6);
        assert_eq!(poly4.linf_norm(), 0);
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    pub fn generate_rq_from_zq_vector(vec: Vec<Zq>) -> Rq {
        let mut temp = [Zq::ZERO; Rq::DEGREE];
        // Process excess terms with sign adjustment
        for i in (0..vec.len()).rev() {
            let m = i / Rq::DEGREE;
            let r = i % Rq::DEGREE;
            let sign = if m % 2 == 0 { 1 } else { -1 };
            if sign == 1 {
                temp[r] += vec[i];
            } else {
                temp[r] -= vec[i];
            }
        }
        Rq::new(temp)
    }

    pub mod helper {
        use super::*;
        pub fn padded(prefix: &[Zq]) -> [Zq; Rq::DEGREE] {
            assert!(
                prefix.len() <= Rq::DEGREE,
                "too many coefficients for degree {}",
                Rq::DEGREE
            );

            let mut out = [Zq::ZERO; Rq::DEGREE];
            out[..prefix.len()].copy_from_slice(prefix);
            out
        }

        pub fn rq_from(prefix: &[Zq]) -> Rq {
            Rq {
                coeffs: padded(prefix),
            }
        }
    }

    // Test new() and polynomial creation
    #[test]
    fn test_new_and_create_poly() {
        let coeffs = [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)];
        let poly = helper::rq_from(&coeffs);
        assert_eq!(poly.coeffs, helper::padded(&coeffs));

        // Direct conversion
        let coeffs2 = [Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)];
        let poly_from_vec_direct: Rq = generate_rq_from_zq_vector(coeffs.to_vec());
        assert_eq!(poly_from_vec_direct.coeffs, helper::padded(&coeffs2));
    }

    #[test]
    fn test_into_coeffs_returns_correct() {
        let mut coeffs = [Zq::ZERO; Rq::DEGREE];
        coeffs[0] = Zq::ONE;
        coeffs[0] = Zq::new(20);
        coeffs[0] = Zq::new(3121);
        coeffs[0] = Zq::new(40);
        let poly: Rq = generate_rq_from_zq_vector(coeffs.to_vec());

        let result = Rq::into_coeffs(poly); // or whatever constructor you have
        assert_eq!(result, coeffs);
    }

    // Test addition of polynomials
    #[test]
    fn test_add() {
        // Within bounds
        let poly1: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        let poly2: Rq =
            generate_rq_from_zq_vector(vec![Zq::new(4), Zq::new(3), Zq::new(2), Zq::ONE]);
        let result = &poly1 + &poly2;
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(5), Zq::new(5), Zq::new(5), Zq::new(5)])
        );

        // Outside of bounds
        let poly3: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        let poly4: Rq =
            generate_rq_from_zq_vector(vec![Zq::NEG_ONE, Zq::new(3), Zq::NEG_ONE, Zq::ONE]);
        let result2 = &poly3 + &poly4;
        assert_eq!(
            result2.coeffs,
            helper::padded(&[Zq::ZERO, Zq::new(5), Zq::new(2), Zq::new(5)])
        );
        // Addition with zero polynomial
        let poly5: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        let poly6: Rq = generate_rq_from_zq_vector(vec![Zq::ZERO]);
        let result3 = &poly5 + &poly6;
        assert_eq!(
            result3.coeffs,
            helper::padded(&[Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)])
        );
        // Addition with high coefficients
        let poly7: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::NEG_ONE]);
        let poly8: Rq =
            generate_rq_from_zq_vector(vec![Zq::NEG_ONE, Zq::NEG_ONE, Zq::NEG_ONE, Zq::NEG_ONE]);
        let result3 = &poly7 + &poly8;
        assert_eq!(
            result3.coeffs,
            helper::padded(&[Zq::ZERO, Zq::ONE, Zq::new(2), -Zq::new(2)])
        );
    }

    // Test multiplication of polynomials
    #[test]
    fn test_mul() {
        // Multiplication with wrapping
        let poly1: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::new(2)]);
        let poly2: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE]);
        let result = &poly1 * &poly2;
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(1), Zq::new(2), Zq::new(3), Zq::new(2)])
        );

        // Multiplication with zero polynomial
        let poly3: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::new(2)]);
        let poly4: Rq = generate_rq_from_zq_vector(vec![Zq::ZERO]);
        let result2 = &poly3 * &poly4;
        assert_eq!(
            result2.coeffs,
            helper::padded(&[Zq::ZERO, Zq::ZERO, Zq::ZERO])
        );

        // Needs to be revised later
        // // Multiplication with wrapping higher order
        // let poly5: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::new(2)]);
        // let poly6: Rq = vec![Zq::ONE, Zq::ONE, Zq::new(7), Zq::new(5)].into();
        // let result3 = poly5 * poly6;
        // assert_eq!(
        //     result3.coeffs,
        //     helper::padded(&[Zq::new(u32::MAX - 12), Zq::new(u32::MAX - 16), Zq::ZERO])
        // );
    }

    // Test subtraction of polynomials
    #[test]
    fn test_sub() {
        // within bounds
        let poly1: Rq =
            generate_rq_from_zq_vector(vec![Zq::new(5), Zq::new(10), Zq::new(15), Zq::new(20)]);
        let poly2: Rq =
            generate_rq_from_zq_vector(vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)]);
        let result = &poly1 - &poly2;
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(3), Zq::new(6), Zq::new(9), Zq::new(12)])
        );

        // Outside of bounds
        let poly3: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ONE, Zq::new(3), Zq::new(2)]);
        let poly4: Rq =
            generate_rq_from_zq_vector(vec![Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)]);
        let result2 = &poly3 - &poly4;
        assert_eq!(
            result2.coeffs,
            helper::padded(&[
                Zq::ZERO - Zq::new(1),
                Zq::ZERO - Zq::new(3),
                Zq::ZERO - Zq::new(3),
                Zq::ZERO - Zq::new(6),
            ])
        );
        // Subtraction with zero polynomial
        let poly5: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        let poly6: Rq = generate_rq_from_zq_vector(vec![Zq::ZERO]);
        let result3 = &poly6 - &poly5;
        let result4 = &poly5 - &poly6;
        assert_eq!(
            result3.coeffs,
            helper::padded(&[
                Zq::ZERO - Zq::new(1),
                Zq::ZERO - Zq::new(2),
                Zq::ZERO - Zq::new(3),
                Zq::ZERO - Zq::new(4),
            ])
        );
        assert_eq!(
            result4.coeffs,
            helper::padded(&[Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)])
        );
    }

    // Test scalar multiplication
    #[test]
    fn test_scalar_mul() {
        let poly: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        let result = &poly * &Zq::new(2);
        assert_eq!(
            result.coeffs,
            helper::padded(&[Zq::new(2), Zq::new(4), Zq::new(6), Zq::new(8)])
        );
    }

    // Test coefficient extraction
    #[test]
    fn test_get_coefficient() {
        let poly: Rq = generate_rq_from_zq_vector(vec![Zq::ONE, Zq::ZERO, Zq::new(5), Zq::NEG_ONE]);
        let vec = helper::padded(&[Zq::ONE, Zq::ZERO, Zq::new(5), Zq::NEG_ONE]).to_vec();
        assert!(poly.coeffs().to_vec() == vec);

        let poly_zero: Rq =
            generate_rq_from_zq_vector(vec![Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO]);
        let vec_zero = helper::padded(&[Zq::ZERO, Zq::ZERO, Zq::ZERO, Zq::ZERO]).to_vec();
        assert!(poly_zero.coeffs().to_vec() == vec_zero);
    }

    #[test]
    fn test_conjugate_automorphism() {
        use crate::core::inner_product::compute_linear_combination;

        let poly1 = helper::rq_from(&[Zq::ONE, Zq::TWO, Zq::new(3)]);
        let poly2 = helper::rq_from(&[Zq::new(4), Zq::new(5), Zq::new(6)]);
        let inner_12 = compute_linear_combination(poly1.coeffs(), poly2.coeffs());
        let conjugated_1 = poly1.conjugate_automorphism();
        let inner_conjugated_12 = &conjugated_1 * &poly2;

        assert_eq!(inner_conjugated_12.coeffs.len(), Rq::DEGREE);
        assert_eq!(inner_conjugated_12.coeffs()[0], Zq::new(32));
        assert_eq!(inner_conjugated_12.coeffs()[1], Zq::new(17));
        assert_eq!(inner_conjugated_12.coeffs()[2], Zq::new(6));

        // ct<\sigma_{-1}(poly1), poly2> ?= <poly1, poly2>
        let ct_inner_conjugated_12 = inner_conjugated_12.coeffs()[0];
        assert_eq!(ct_inner_conjugated_12, inner_12);
    }
}

#[cfg(test)]
mod decomposition_tests {
    use crate::ring::rq::tests::{generate_rq_from_zq_vector, helper};

    use super::*;

    #[test]
    fn failure_test() {
        let poly1 =
            generate_rq_from_zq_vector(vec![-Zq::new(192), -Zq::new(0), -Zq::new(0), -Zq::new(0)]);
        let decomposed = poly1.decompose(Zq::new(1802), 10);

        let mut result = Rq::zero();
        let mut exponential_base = Zq::new(1);
        for elem in decomposed {
            result = &result + &(&elem * &exponential_base);
            exponential_base *= Zq::new(24);
        }
        assert_eq!(result, poly1);
    }

    #[test]
    fn test_base2_decomposition() {
        // Test case 1: Base 2 decomposition
        let poly: Rq =
            generate_rq_from_zq_vector(vec![Zq::new(5), Zq::new(3), Zq::new(7), Zq::new(1)]);
        let parts = poly.decompose(Zq::TWO, 4);

        // Part 0: remainders mod 2 (no centering needed for base 2)
        assert_eq!(
            parts[0].coeffs,
            helper::padded(&[
                Zq::ONE, // 5 mod 2 = 1
                Zq::ONE, // 3 mod 2 = 1
                Zq::ONE, // 7 mod 2 = 1
                Zq::ONE, // 1 mod 2 = 1
            ])
        );

        // Part 1: quotients after division by 2
        assert_eq!(
            parts[1].coeffs,
            helper::padded(&[Zq::ZERO, Zq::ONE, Zq::ONE, Zq::ZERO,])
        );

        // Part 2
        assert_eq!(
            parts[2].coeffs,
            helper::padded(&[
                Zq::ONE,  // 5 div 2 = 2
                Zq::ZERO, // 3 div 2 = 1
                Zq::ONE,  // 7 div 2 = 3
                Zq::ZERO, // 1 div 2 = 0
            ])
        );

        // Verify Base 2 reconstruction coefficient by coefficient
        let mut result = Rq::zero();
        let mut expo_base = Zq::ONE;
        for part in parts {
            result = &result + &(&part * &expo_base);
            expo_base *= Zq::TWO;
        }
        assert_eq!(poly, result);
    }

    #[test]
    fn test_decomposition_edge_cases() {
        // Test zero polynomial
        let zero_poly: Rq = generate_rq_from_zq_vector(vec![Zq::ZERO; 4]);
        let parts = zero_poly.decompose(Zq::TWO, 2);
        assert!(
            // Check any polynomial is zero
            parts
                .iter()
                .all(|p| p.coeffs().iter().all(|&coeff| coeff == Zq::ZERO)),
            "Zero polynomial decomposition failed"
        );

        // Test single part decomposition
        let simple_poly: Rq =
            generate_rq_from_zq_vector(vec![Zq::ONE, Zq::new(2), Zq::new(3), Zq::new(4)]);
        let parts = simple_poly.decompose(Zq::TWO, 1);
        assert_eq!(parts.len(), 1, "Single part decomposition length incorrect");
    }

    #[test]
    fn test_decomposition_real_values() {
        let poly1 = Rq {
            coeffs: [
                Zq::new(23071),
                Zq::new(4294941461),
                Zq::new(4084),
                Zq::new(23413),
                Zq::new(26902),
                Zq::new(4294960817),
                Zq::new(4294957829),
                Zq::new(18248),
                Zq::new(15293),
                Zq::new(45978),
                Zq::new(27197),
                Zq::new(11233),
                Zq::new(4294962536),
                Zq::new(1196),
                Zq::new(17083),
                Zq::new(4294960822),
                Zq::new(4294967103),
                Zq::new(4294949429),
                Zq::new(2926),
                Zq::new(2742),
                Zq::new(27552),
                Zq::new(4294955871),
                Zq::new(2917),
                Zq::new(4294938904),
                Zq::new(2288),
                Zq::new(4294948781),
                Zq::new(4294966588),
                Zq::new(4294951307),
                Zq::new(4294961210),
                Zq::new(19355),
                Zq::new(4294956913),
                Zq::new(14940),
                Zq::new(4294934068),
                Zq::new(4294961161),
                Zq::new(4294959017),
                Zq::new(3339),
                Zq::new(4294932967),
                Zq::new(4294938251),
                Zq::new(4294943006),
                Zq::new(4294965167),
                Zq::new(4294960800),
                Zq::new(4161),
                Zq::new(9213),
                Zq::new(4294962556),
                Zq::new(14299),
                Zq::new(44418),
                Zq::new(2438),
                Zq::new(6583),
                Zq::new(4294957783),
                Zq::new(16),
                Zq::new(4294941724),
                Zq::new(4294966309),
                Zq::new(4294966984),
                Zq::new(4294956138),
                Zq::new(1779),
                Zq::new(29598),
                Zq::new(16393),
                Zq::new(728),
                Zq::new(4294944371),
                Zq::new(4294951792),
                Zq::new(4294943824),
                Zq::new(4294937618),
                Zq::new(4294955208),
                Zq::new(11235),
            ],
        };
        let base = Zq::new(1802);
        let decomposed = poly1.decompose(base, 2);

        // Verify reconstruction
        let reconstructed = &decomposed[0] + &(&decomposed[1] * &base);

        assert_eq!(decomposed.len(), 2);
        assert_eq!(reconstructed, poly1, "Base {}: Reconstruction failed", base);
    }

    #[test]
    fn test_linf_norm() {
        let poly1 = Rq {
            coeffs: [
                Zq::new(23071),
                Zq::new(4294941461),
                Zq::new(4084),
                Zq::new(23413),
                Zq::new(26902),
                Zq::new(4294960817),
                Zq::new(4294957829),
                Zq::new(18248),
                Zq::new(15293),
                Zq::new(45978),
                Zq::new(27197),
                Zq::new(11233),
                Zq::new(4294962536),
                Zq::new(1196),
                Zq::new(17083),
                Zq::new(4294960822),
                Zq::new(4294967103),
                Zq::new(4294949429),
                Zq::new(2926),
                Zq::new(2742),
                Zq::new(27552),
                Zq::new(4294955871),
                Zq::new(2917),
                Zq::new(4294938904),
                Zq::new(2288),
                Zq::new(4294948781),
                Zq::new(4294966588),
                Zq::new(4294951307),
                Zq::new(4294961210),
                Zq::new(19355),
                Zq::new(4294956913),
                Zq::new(14940),
                Zq::new(4294934068),
                Zq::new(4294961161),
                Zq::new(4294959017),
                Zq::new(3339),
                Zq::new(4294932967),
                Zq::new(4294938251),
                Zq::new(4294943006),
                Zq::new(4294965167),
                Zq::new(4294960800),
                Zq::new(4161),
                Zq::new(9213),
                Zq::new(4294962556),
                Zq::new(14299),
                Zq::new(44418),
                Zq::new(2438),
                Zq::new(6583),
                Zq::new(4294957783),
                Zq::new(16),
                Zq::new(4294941724),
                Zq::new(4294966309),
                Zq::new(4294966984),
                Zq::new(4294956138),
                Zq::new(1779),
                Zq::new(29598),
                Zq::new(16393),
                Zq::new(728),
                Zq::new(4294944371),
                Zq::new(4294951792),
                Zq::new(4294943824),
                Zq::new(4294937618),
                Zq::new(4294955208),
                Zq::new(11235),
            ],
        };
        let base = Zq::new(1802);
        let decomposed = poly1.decompose(base, 2);

        // Verify reconstruction
        for decomposed_poly in decomposed {
            assert!(decomposed_poly.linf_norm() <= 901)
        }
    }
}
