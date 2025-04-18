# Changing the Modulus & Norm Checks
This document clarifies handling different moduli and verifying norm constraints in the process of adapting Falcon signatures for aggregation within the LaBRADOR proof system.


## Different Moduli
Falcon and LaBRADOR both operate over polynomial rings of the form $\mathcal{R}\_q = \mathbb{Z}\_q[X]/(X^d + 1)$, but they use different moduli: Falcon uses a smaller modulus $q$, while LaBRADOR requires a larger modulus $q' > q$. This difference arises because Falcon’s standard parameter sets do not meet the conditions for the Johnson-Lindenstrauss projection used in LaBRADOR. Therefore, we must lift Falcon signatures verification into a larger ring modulo $q'$, while ensuring compatibility and security.


## Falcon Signature Scheme
In Falcon, a signer uses their secret key to generate a signature consisting of two small polynomials $\mathbf{s}\_{i\_1}, \mathbf{s}\_{i\_2} \in \mathcal{R}\_q$ for a message $m\_i$. For verification, the signature must satisfy two conditions:
1. Verification Equation:

$$\mathbf{s}\_{i\_1} + \mathbf{h}\_i \mathbf{s}\_{i\_2} = \mathbf{t}\_i \mod q \qquad (1)$$
  
Here, $\mathbf{h}\_i$ is the public key (a polynomial in $\mathcal{R}\_q$), $\mathbf{t}\_i = H(r\_i, m\_i)$ is the hash of the message $m\_i$ and a random salt $r\_i \in \\{0, 1\\}^{320}$, and $H$ is a cryptographic hash function.

2. Norm Bound:
   
$$\\|(\mathbf{s}\_{i\_1}, \mathbf{s}\_{i\_2}) \\|\_2 = \sqrt{\\| \mathbf{s}\_{i\_1} \|\_2^2 + \\| \mathbf{s}\_{i\_2} \\|\_2^2} \leq \beta \ll q \qquad (2)$$

The $\ell\_2$-norm bound $\beta$ ensures the signature is small, a critical property for Falcon’s security.

Note that (1) and (2) hold modulo $q$, but not necessarily modulo $q'$ or over the integers ($\mathcal{R}$).


## Falcon Verification in LaBRADOR
To aggregate Falcon signatures in LaBRADOR, we need to adapt the verification process to work over $\mathcal{R}\_{q'} = \mathbb{Z}\_{q'}[X]/(X^d + 1)$. Since equation (1) is defined modulo $q$, we must rewrite it to be compatible with $q'$ while ensuring no wrap-around occurs.

### Rewriting Equation (1) for Modulo $q'$
Equation (1) does not naturally hold over the integers $( \mathcal{R} )$ or modulo $q'$. To address this, we introduce an additional witness polynomial $\mathbf{v}\_i \in \mathcal{R}\_{q'}$ and reformulate (1) as:

$$\mathbf{s}\_{i\_1} + \mathbf{h}\_i \mathbf{s}\_{i\_2} + q \mathbf{v}\_i - \mathbf{t}\_i = \mathbf{0} \in \mathcal{R} \qquad (3) $$

This equation holds over the integers by construction, as $q \mathbf{v}\_i$ accounts for the difference modulo $q$. However, LaBRADOR operates over $\mathcal{R}\_{q'}$, so we need:

$$\mathbf{s}\_{i\_1} + \mathbf{h}\_i \mathbf{s}\_{i\_2} + q \mathbf{v}\_i - \mathbf{t}\_i = \mathbf{0} \in \mathcal{R}\_{q'} \qquad (4)$$

For (4) to be equivalent to (3), the coefficients of the left-hand side must not wrap around modulo $q'$. This requires the infinity norm of the sum to be less than $q'/2$:

$$\\| \mathbf{s}\_{i\_1} + \mathbf{h}\_i \mathbf{s}\_{i\_2} + q \mathbf{v}\_i - \mathbf{t}\_i \\|\_\infty < \frac{q'}{2} $$

To ensure this, we impose bounds on each term:

$$\\| \mathbf{s}\_{i\_1} \\|\_\infty < \frac{q'}{6} \qquad 
\\| \mathbf{h}\_i \mathbf{s}\_{i\_2} \\|\_\infty < \frac{q'}{6} \qquad
\\| q \mathbf{v}\_i \\|\_\infty < \frac{q'}{6}$$

This implies

$$\\| \mathbf{s}\_{i\_1} \\|\_\infty < \frac{q'}{6} \qquad 
\\| \mathbf{s}\_{i\_2} \\|\_\infty < \frac{q'}{6qd} \qquad
\\| q \mathbf{v}\_i \\|\_\infty < \frac{q'}{6q}$$

These bounds ensure the total sum stays within $\frac{q'}{2}$, preventing wrap-around. In LaBRADOR, $\mathbf{s}\_{i\_1}, \mathbf{s}\_{i\_2}, \mathbf{v}\_i$ are witnesses, and we prove (4) holds.


## Verifying the Norm Bound in LaBRADOR
Falcon’s norm condition (2) must also be proven in LaBRADOR without revealing $\mathbf{s}\_{i\_1}$ or $\mathbf{s}\_{i\_2}$. We rewrite:

$$ \beta^2 - \\| \mathbf{s}\_{i\_1} \\|\_2^2 - \\| \mathbf{s}\_{i\_2} \\|\_2^2 \geq 0 $$

By Lagrange’s four-square theorem, any non-negative integer is the sum of four squares. Thus, we introduce four integer witnesses $\epsilon\_{i,0}, \epsilon\_{i,1}, \epsilon\_{i,2}, \epsilon\_{i,3}$ such that:

$$ \beta^2 - \\| \mathbf{s}\_{i\_1} \\|\_2^2 - \\| \mathbf{s}\_{i\_2} \\|\_2^2 = \epsilon\_{i,0}^2 + \epsilon\_{i,1}^2 + \epsilon\_{i,2}^2 + \epsilon\_{i,3}^2 \qquad (5)$$

To fit LaBRADOR’s polynomial constraints, we represent each $\epsilon\_{i,j}$ as a polynomial of degree at most 3

$$ \mathbf{\epsilon}\_i = \epsilon\_{i,0} + \epsilon\_{i,1} X + \epsilon\_{i,2} X^2 + \epsilon\_{i,3} X^3 $$

The norm squared is computed using the conjugation automorphism $\sigma\_{-1}$, where $\\| \mathbf{v} \\|\_2^2 = \sigma\_{-1}(\mathbf{v}) \mathbf{v}$. We rewrite (5) as:

$$ \sigma\_{-1}(\mathbf{\epsilon}\_i) \mathbf{\epsilon}\_i = \beta^2 - \sigma\_{-1}(\mathbf{s}\_{i\_1}) \mathbf{s}\_{i\_1} - \sigma\_{-1}(\mathbf{s}\_{i\_2}) \mathbf{s}\_{i\_2} \in \mathcal{R} $$

In LaBRADOR, this must hold modulo $q'$:

$$ \sigma\_{-1}(\mathbf{\epsilon}\_i) \mathbf{\epsilon}\_i - (\beta^2 - \sigma\_{-1}(\mathbf{s}\_{i\_1}) \mathbf{s}\_{i\_1} - \sigma\_{-1}(\mathbf{s}\_{i\_2}) \mathbf{s}\_{i\_2}) = 0 \mod q' \qquad (6) $$

To avoid wrap-around, the infinity norm of the left-hand side must be less than $q'/2$:

$$ \\| \sigma\_{-1}(\mathbf{\epsilon}\_i) \mathbf{\epsilon}\_i - (\beta^2 - \sigma\_{-1}(\mathbf{s}\_{i\_1}) \mathbf{s}\_{i\_1} - \sigma\_{-1}(\mathbf{s}\_{i\_2}) \mathbf{s}\_{i\_2}) \\|\_\infty < \frac{q'}{2} $$

This involves $2d + 4$ coefficients ( $2d$ from $\mathbf{s}\_{i\_1}$, $\mathbf{s}\_{i\_2}$, and 4 from $\mathbf{\epsilon}\_i $). Bounding each coefficient’s $\ell\_\infty$-norm by $\sqrt{\frac{q'}{2(2d + 4)}}$ ensures the sum stays within $q'/2$.


### Conjugation Automorphism Checks

For a polynomial $\mathbf{f}(X) = a\_0 + a\_1 X + \cdots + a\_{d-1} X^{d-1}$, the conjugation automorphism is $\sigma\_{-1}(\mathbf{f}) = a\_0 + a\_1 X^{-1} + \cdots + a\_{d-1} X^{-(d-1)}$. In $\mathcal{R}\_{q'}$, $X^{-t}$ is computed modulo $X^d + 1$. We verify coefficients using dot product constraints:
1. Check the $j$-th coefficient equals a value $b$:

$$ \sigma\_{-1}(X^j) \mathbf{a} = b \mod q'$$


2. Check two coefficients $a\_j, c\_k$ match:

$$ \sigma\_{-1}(X^j) \mathbf{a} = \sigma\_{-1}(X^k) \mathbf{c} \mod q'$$

For $\mathbf{\epsilon}\_i$, we ensure the degree is at most 3 by setting coefficients $\epsilon\_{i,4}, \ldots, \epsilon\_{i,d-1}$ to zero. For $\mathbf{s}\_{i\_1}, \mathbf{s}\_{i\_2}, \mathbf{\epsilon}\_i$, we verify the conjugation and coefficients match in (6). This requires approximately $4d n$ extra constraints for $n$ signatures.



## Choosing the Modulus $q'$
To satisfy all constraints, $q'$ must be sufficiently large. After computations for soundness and completeness:

$$ q' > \frac{1024}{15} (d + 2) \beta^2 N$$

- Falcon-512: $q' > 2^{40.12} N$
- Falcon-1024: $q' > 2^{42.16} N$

For $N = 2^{20}$ signatures:
- Falcon-512: $q' > 2^{60.12}$, so a 61-bit modulus suffices.
- Falcon-1024: $q' > 2^{62.16}$, so a 63-bit modulus is needed.


## Conclusion
This guide clarifies how to adapt Falcon signatures for LaBRADOR by changing the modulus and verifying norms. By rewriting the verification equation, introducing witnesses, and enforcing norm bounds, developers can implement a secure and correct aggregation system.