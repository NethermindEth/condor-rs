# Changing the Modulus & Norm Checks

Although Falcon and LaBRADOR were designed over a similar structure of a polynomial ring $\mathcal{R}\_{q} = Z_{q}\[ x \]/(x^{d}-1)$ of degree $d$ and modulus $q$, they both use different moduli. We will denote $q$ as the Falcon modulus and $q' > q$ as the LaBRADOR one. When using Falcon to sign a message $m$, the signer uses its secret key to obtain two small lattice vectors $(\mathbf{s}\_{i1}, \mathbf{s}\_{i2})$ such that:
$$\mathbf{s}\_{i1}+\mathbf{hs}\_{i2} = H(r,m) \mod q$$
$$\lVert (\mathbf{s}\_{i1}, \mathbf{s}\_{i2}) \rVert_{2} \leq \beta$$
Where $H()$ is a hash function, $\mathbf{h}$ is the public key, and $r \in \\{0, 1\\}^{320}$ is a random salt.

It's important to notice these equations are valid mod $q$, which means they may not be valid mod $q'$. So, in order to make an equivalent version of these constraints but in the LaBRADOR ring, we will find an equivalent restriction over $\mathcal{R}$. Since the validity of a restriction being equal to zero in $\mathcal{R}$ should still hold under some modulus wrapping in $\mathcal{R}\_{q^{'}}$. In order for the restrictions to hold in $\mathcal{R}$, we would like intuitively that $(\mathbf{s}\_{i1}, \mathbf{s}\_{i2})$ will be "small enough" that they would not wrap around over the $q'$ modulus. Since we are not in control of the size of $(\mathbf{s}\_{i1}, \mathbf{s}\_{i2})$, the main idea behind this will be to slightly modify the restrictions from Falcon so that we can guarantee that for a large enough $q'$, they will still be valid over this new ring.

For the first restriction, we know from Falcon that it must be valid in $\mathcal{R}\_{q}$. This isn't necessarily valid in $\mathcal{R}$. In order to force it to be valid, we will add an extra witness $v_{i} \in \mathcal{R}\_{q'}$, leaving us with this restriction over $\mathcal{R}$:
$$\mathbf{s}\_{i1}+\mathbf{hs}\_{i2}+qv_{i} - H(r,m) = 0$$

Since adding multiples of $q$ to the restriction doesn't affect the original formulation $\mod q$, then this restriction continues to be valid in $\mathcal{R}\_{q}$. In order to answer how big $q'$ needs to be so that the witness vectors (including the new one) don't wrap around, we need to prove that the witnesses are small enough.

Now, in order to find an equivalent restriction for the restriction on the size of the witnesses, we will rewrite $\lVert (\mathbf{s}\_{i1}, \mathbf{s}\_{i2}) \rVert_{2} \leq \beta$ as:
$$\lVert \mathbf{s}\_{i1}\rVert^{2}+ \lVert \mathbf{s}\_{2}\rVert^{2} \leq \beta^{2}$$
Where $\lVert . \rVert$ is still the Euclidean norm. Because we know that $\beta^{2} - \lVert \mathbf{s}\_{i1}\rVert^{2} - \lVert \mathbf{s}\_{i2}\rVert^{2}$ is non-negative, we can use Lagrange’s four-square theorem to rewrite it as a sum of four squared numbers $\epsilon_{0}^{2}, \epsilon_{1}^{2}, \epsilon_{2}^{2}, \epsilon_{3}^{2}$.

Because all dot product constraints in LaBRADOR are of the form:
$$
f(\mathbf{\bar{\mathbf{s}}})=0 \text{ or } ct(f(\mathbf{\bar{\mathbf{s}}})) = 0
$$
we can rewrite now the restriction on the size of the witnesses in the accepted format by writing $\epsilon_{0}, \epsilon_{1}X, \epsilon_{2}X^{2}, \epsilon_{3}X^{3} = \epsilon$ and:
$$ ct(\sigma_{-1}(\epsilon)\epsilon - (\beta^{2}-\sigma_{-1}(\mathbf{s}\_{i1})\mathbf{s}\_{i1}-\sigma_{-1}(\mathbf{s}\_{i2})\mathbf{s}\_{i2})) = 0 \mod q' $$
where $\sigma_{-1}$ is the conjugation automorphism. The validity of this dot constraint over $q'$ tells us the bound is valid since the existence of the polynomial $\epsilon$ that follows such constraints implies that $\beta^{2} - \lVert \mathbf{s}\_{i1}\rVert^{2} - \lVert \mathbf{s}\_{i2}\rVert^{2}$ is non-negative.

Besides this new dot constraint, in order to show that both the automorphisms and the epsilon polynomial were correctly calculated inside the constraint, we will need to add small extra dot product constraints.

For the conjugated automorphism, since it involves a permutation of the polynomial's coefficients and a change of sign, it will suffice to check for each element that the permutation was done correctly as well as the sign change. In the case of the epsilon polynomial, one would only need to check the degree is at most $4$ by checking all other coefficients $\epsilon_{4} \dots \epsilon_{d-1}$ to be zero. For both cases, checking an element of a polynomial is equal to some element $b$ or to some element from another polynomial $c$ can be written as a dot product in this way: 

$$ct(\sigma_{-1}(X^{j})a -b) = 0 \mod q'$$
$$ct(\sigma_{-1}(X^{j})a -\sigma_{-1}(X^{k})c) = 0 \mod q'$$




 

