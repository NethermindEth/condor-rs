//! Provides utility for decomposition parameters. Adapted from stark-rings:
//! https://github.com/NethermindEth/stark-rings/tree/main/ring/src/balanced_decomposition
//!
//! Section 5.2 Main Protocol, decompose \vec{t_i} into t_1 >= 2 parts
//! with a small base b_1, where each coefficient's infity norm should
//! be bounded to b_1/2.
use crate::{rq::Rq, rq_vector::RqVector, zq::Zq};
pub trait Decompose: Sized + Copy {
    fn decompose_balanced_in_place(&self, b: Zq, out: &mut [Zq]);
    fn decompose(&self, b: Zq, padding_size: usize) -> Vec<Zq> {
        let mut res = vec![Zq::ZERO; padding_size];

        self.decompose_balanced_in_place(b, &mut res);

        res
    }
}

/// Returns the balanced decomposition of a slice as a Vec of Vecs,
/// "Balanced" means that each coefficient will lie in the range roughly
/// [-b/2, b/2], from the very top of Page 12
///
/// # Arguments
/// * `v`: input coefficient
/// * `b`: basis for the decomposition
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`,
///   or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
pub fn decompose_balanced_in_place(v: Zq, b: Zq, out: &mut [Zq]) {
    assert!(
        b != Zq::ONE && b != Zq::ZERO,
        "cannot decompose in basis 0 or 1"
    );

    let mut curr = v.to_signed_zq();
    let mut current_i = 0;
    let b = b.to_i128();
    let b_half_floor = b.div_euclid(2);
    loop {
        let rem = curr % b; // rem = curr % b is in [-(b-1), (b-1)]

        // Ensure digit is in [-b/2, b/2]
        if rem.abs() <= b_half_floor {
            out[current_i] = Zq::to_unsigned_zq(rem);
            curr /= b; // Rust integer division rounds towards zero
        } else {
            // The next element in the decomposition is sign(rem) * (|rem| - b)
            if rem < 0.into() {
                out[current_i] = Zq::to_unsigned_zq(rem + b);
            } else {
                out[current_i] = Zq::to_unsigned_zq(rem - b);
            }
            let carry = rounded_div(rem, b); // Round toward nearest integer, not towards 0
            curr = (curr / b) + carry;
        }

        current_i += 1;

        if curr == 0 {
            break;
        }
    }

    for out_tail_elem in out[current_i..].iter_mut() {
        *out_tail_elem = Zq::ZERO;
    }
}

/// Divides `n` by `d`, rounding to the nearest integer.
/// If the fractional part is exactly 0.5, the function rounds away from zero.
///
/// # Panics
///
/// Panics if `d` is zero.
fn rounded_div(n: i128, d: i128) -> i128 {
    assert!(d != 0, "division by zero");

    let q = n / d; // Truncated quotient (rounds toward zero)
    let r = n % d; // Remainder

    // If the absolute value of the remainder is at least half of |d|,
    // adjust the quotient in the direction of n's sign.
    if 2 * r.abs() >= d.abs() {
        if n >= 0 {
            q + 1
        } else {
            q - 1
        }
    } else {
        q
    }
}

/// Returns the balanced decomposition of a slice as a Vec of Vecs.
pub fn decompose_balanced_vec<const D: usize>(
    v: Rq<D>,
    b: Zq,
    padding_size: usize,
) -> Vec<Vec<Zq>> {
    v.iter()
        .map(|v_i| decompose_balanced(*v_i, b, padding_size))
        .collect()
}

pub fn decompose_balanced_commits<const N: usize, const D: usize>(
    v: RqVector<N, D>,
    b: Zq,
    padding_size: usize,
) -> Vec<Vec<Vec<Zq>>> {
    v.as_slice()
        .iter()
        .map(|v_i| {
            v_i.iter()
                .map(|v_j| decompose_balanced(*v_j, b, padding_size))
                .collect()
        })
        .collect()
}

pub fn decompose_balanced(v: Zq, b: Zq, padding_size: usize) -> Vec<Zq> {
    let mut result = vec![Zq::ZERO; padding_size];

    decompose_balanced_in_place(v, b, &mut result);

    result
}

/// revert the decomposition, mainly for testing
pub fn recompose(decomp: Vec<Zq>, b: Zq) -> Zq {
    let mut result = Zq::ZERO;

    for v_i in decomp.iter().rev() {
        result *= b;
        result += *v_i;
    }

    result
}

#[cfg(test)]
pub mod test_coeffs {
    use crate::rq::Rq;
    use crate::zq::Zq;

    use super::*;

    const BASIS_TEST_RANGE: [Zq; 1] = [Zq::new(4)];

    #[test]
    fn test_decompose_balanced() {
        let poly1: Rq<2> = vec![Zq::new(100), Zq::new(200)].into();
        const D: usize = 8;
        let b = Zq::new(2);
        let b_half = b / Zq::new(2);
        let elements = poly1.get_coefficients();
        let decomp1: Vec<Zq> = decompose_balanced(elements[0], b, D);
        let decomp2: Vec<Zq> = decompose_balanced(elements[1], b, D);

        for v1_i in &decomp1 {
            assert!(*v1_i <= b_half || *v1_i >= -b_half);
        }

        for v2_i in &decomp2 {
            assert!(*v2_i <= b_half || *v2_i >= -b_half);
        }
        assert_eq!(elements[0], recompose(decomp1, b));
        assert_eq!(elements[1], recompose(decomp2, b));
    }

    #[test]
    fn test_decompose_balanced_vec() {
        let poly1: Rq<2> = vec![Zq::new(100), Zq::new(200)].into();
        for b in BASIS_TEST_RANGE {
            let b_half = b / Zq::new(2);
            let decomp = decompose_balanced_vec(poly1.clone(), b, 8);
            println!("{:?}", decomp);

            // Check that all entries are smaller than b/2 in absolute value
            for d_i in &decomp {
                for d_ij in d_i {
                    assert!(*d_ij <= b_half || *d_ij >= -b_half);
                }
            }
            println!("{:?}", decomp.len());
            for (i, item) in decomp.iter().enumerate() {
                // Check that the decomposition is correct
                let decomp_i = item.clone();
                println!("{:?}", decomp_i);
                assert_eq!(poly1.get_coefficients()[i], recompose(decomp_i, b));
            }
        }
    }

    #[test]
    fn test_decompose_balanced_commits() {
        let poly5: Rq<2> = vec![Zq::new(100), Zq::new(200)].into();
        let poly6: Rq<2> = vec![Zq::new(1000), Zq::new(2000)].into();
        // let poly7: Rq<2> = vec![Zq::new(10000), Zq::new(20000)].into();
        // let poly8: Rq<2> = vec![Zq::new(100000), Zq::new(200000)].into();
        let vec_5: RqVector<2, 2> = RqVector::from(vec![poly5, poly6]);
        // let vec_6: RqVector<2, 2> = RqVector::from(vec![poly7, poly8]);
        for b in BASIS_TEST_RANGE {
            let b_half = b / Zq::new(2);
            let decomp = decompose_balanced_commits(vec_5.clone(), b, 8);
            println!("{:?}", decomp);

            // Check that all entries are smaller than b/2 in absolute value
            for d_i in &decomp {
                for d_ij in d_i {
                    for d_ijk in d_ij {
                        assert!(*d_ijk <= b_half || *d_ijk >= -b_half);
                    }
                }
            }
            for i in 0..decomp.len() {
                // Check that the decomposition is correct
                let decomp_i = decomp[i].clone();
                for (j, item) in decomp_i.iter().enumerate() {
                    let decomp_ij = item.clone();
                    assert_eq!(vec_5[i].get_coefficients()[j], recompose(decomp_ij, b));
                }
            }
        }
    }
}
