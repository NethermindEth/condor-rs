#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Poly<const D: usize> {
    coeffs: [u32; D],
}

impl<const D: usize> Poly<D> {
    // Constructor for the polynomial ring
    pub fn new(coeffs: [u32; D]) -> Self {
        Poly { coeffs }
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut result = [0u32; D];
        for (r, (a, b)) in result
            .iter_mut()
            .zip(self.coeffs.iter().zip(other.coeffs.iter()))
        {
            *r = a.wrapping_add(*b);
        }
        Poly::new(result)
    }

    // Create a Polynomial from a vector
    pub fn create_poly(coeffs: Vec<i32>) -> Poly<D> {
        let mut arr = [0u32; D];
        let u32_coeffs: Vec<u32> = coeffs.iter().map(|&coeff| coeff as u32).collect();
        // First D elements are assigned directly
        for (i, &u32_coeffs) in u32_coeffs.iter().take(D).enumerate() {
            arr[i] = u32_coeffs;
        }

        // Handle additional elements by subtracting them at (index % D)
        for (i, &u32_coeffs) in u32_coeffs.iter().skip(D).enumerate() {
            let mod_index = i % D;
            arr[mod_index] = arr[mod_index].wrapping_sub(u32_coeffs);
        }

        Poly::new(arr)
    }

    // Polynomial multiplication modulo x^D + 1
    pub fn mul(&self, other: &Self) -> Self {
        let mut result = [0u32; D];

        for i in 0..D {
            for j in 0..D {
                let degree = (i + j) % D;
                if (i + j) > D {
                    result[degree] =
                        result[degree].wrapping_sub(self.coeffs[i].wrapping_mul(other.coeffs[j]));
                } else {
                    // normal multiplication
                    result[degree] =
                        result[degree].wrapping_add(self.coeffs[i].wrapping_mul(other.coeffs[j]));
                }
            }
        }

        Poly::new(result)
    }

    // Dot product between coefficients
    pub fn inner_product(&self, other: &Self) -> u32 {
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .map(|(&a, &b)| a.wrapping_mul(b))
            .fold(0u32, u32::wrapping_add)
    }
}
