# Changing the Modulus & Norm Checks

## Different Modulos
While Falcon and LaBRADOR are designed over a similar structure of a polynomial rings $\mathcal{R}\_{q} = Z_{q}[ x ]/(x^{d}-1)$, they use different modulies. Falcon uses modulo $q$ and LaBRADOR uses modulo $q'$, while $q' > q$. The  reason is that nighther of the parameter sets of Falcon satisfy the condition for Johnson-Lindenstrauss projection of LaBRADOR. So we need to use a larger modulo. 

#### Modulo in Falcon Signature Scheme
To sign a message $m$ in Falcon, the signer uses its secret key to obtain two small polynomials $(\mathbf{s}_{i_1}, \mathbf{s}_{i_2})$ such that:
$$\mathbf{s}_{i_1}+\mathbf{hs}_{i_2} = H(r,m) = \mathbf{t}_i \mod q \qquad (1)$$

$$\qquad \| (\mathbf{s}_{i_1}, \mathbf{s}_{i_2}) \|_{2} \leq \beta, \qquad (2)$$
Where $H()$ is a hash function, $\mathbf{h}$ is the public key, and $r \in \{0, 1\}^{320}$ is a random salt.

It is important to notice these equations are valid mod $q$, and they are not neccesarily  valid mod $q'$. 

#### Falcon Verification in LaBRADOR Modulo
To work with LaBRADOR, we need to work modulo $q'$, so we need to define an equivalent formula of Falcon signauture verification (equation (1)), but in moulo $q'$.

However, equation (1) is not neccessarily valid in $\mathcal{R}$ or $\mathcal{R}_{q'}$. We know that is an equation is valid in $\mathcal{R}$, it is clearly valid in $\mathcal{R}_{q}$ and $\mathcal{R}_{q'}$, but it is not correct for the other side.

The idea to resolve this issue is that, first convert the equation (1) to a valid equation holds in $\mathcal{R}$, then, by forcing small norms, there would be no wrap-around modulo $q'$, and the equation can bw represented in $\mathcal{R}_{q'}$.

First, (1) is not neccesarily valid in $\mathcal{R}$. In order to force it to be valid, we will add an extra witness $v_i \in \mathcal{R}_{q'}$, and rewrite formula (1) as follow, which holds for any ring $\mathcal{R}$:
$$ \mathbf{s}_{i_1}+\mathbf{hs}_{i_2} + q\mathbf{v}_i - \mathbf{t}_i = \mathbf{0} \in \mathcal{R} \qquad (3) $$

LaBRADOR constraints are defined over $\mathcal{R}_{q'}$, and equation (3) is not neccesarily valid in $\mathcal{R}_{q'}$. 
However, if these polynomial functions are so small that could never cause a wrap-around modulo $q'$, we can write (3) as:
$$ \mathbf{s}_{i_1}+\mathbf{hs}_{i_2} + q\mathbf{v}_i - \mathbf{t}_i = \mathbf{0} \in \mathcal{R}_{q'} \qquad (4) $$


### Norm Checks
We need to show that:
1. norm check in (2) is satisfied, 
2. Witnesses $\mathbf{s}_{i_1}, \mathbf{hs}_{i_2}, \mathbf{v}_{i}$ have small norms so that (3) and (4) are equivalent for these norms.

So we proceed as follow:
1.  we can rewrite $\| (\mathbf{s}_{i_1}, \mathbf{s}_{i_2}) \|_{2} \leq \beta$ as:
$$\| \mathbf{s}_{i_1}\|^{2}+ \| \mathbf{s}_{i_2}\|^{2} \leq \beta^{2} \implies \beta^{2} - \| \mathbf{s}_{i_1}\|^{2} - \| \mathbf{s}_{i_2}\|^{2}\ge 0$$
Where $\| . \|$ is still the Euclidean norm.

According to Lagrange's four-square theorem, any non-negative number can be written as sum of four square integers.
Therefore, as $\beta^{2} - \| \mathbf{s}_{i_1}\|^{2} - \| \mathbf{s}_{i_2}\|^{2} \ge 0$ , we can use Lagrange’s four-square theorem to rewrite it as a sum of four squared integers $\epsilon_{0}^{2}, \epsilon_{1}^{2}, \epsilon_{2}^{2}, \epsilon_{3}^{2}$:
$$\beta^{2} - \| \mathbf{s}_{i_1}\|^{2} - \| \mathbf{s}_{i_2}\|^{2} = \epsilon_{0}^{2} + \epsilon_{1}^{2} + \epsilon_{2}^{2} + \epsilon_{3}^{2} \qquad (5)$$

Let $$\epsilon_i = \epsilon_{0}^{2}+ \epsilon_{1}^{2}X + \epsilon_{2}^{2}X^2 + \epsilon_{3}^{2}X^3$$

To make them compatible with LaBRADOR's constraints, which are of the form:
$$
f(\mathbf{\bar{\mathbf{s}}})=0 \text{ or } ct(f(\mathbf{\bar{\mathbf{s}}})) = 0,$$
we can rewrite (5) as:
$$ct(\sigma_{-1}(\epsilon_i)\epsilon_i - (\beta^{2}-\sigma_{-1}(\mathbf{s}_{i_1})\mathbf{s}_{i_1}-\sigma_{-1}(\mathbf{s}_{i_2})\mathbf{s}_{i_2})) = 0 \mod q' \qquad (9)$$
where $\sigma_{-1}$ is the conjugation automorphism. The validity of this dot constraint over $q'$ tells us the bound is valid since the existence of the polynomial $\epsilon$ that follows such constraints implies that $\beta^{2} - \| \mathbf{s}_{i_1}\|^{2} - \| \mathbf{s}_{i_2}\|^{2}$ is non-negative.

Besides this new dot constraint, in order to show that both the automorphisms and the epsilon polynomial were correctly calculated inside the constraint, we will need to add almost $4dn$ extra dot product constraints ($n$ signatures, $d$ coefficients, and 4 functions $\epsilon_i, \sigma_{-1}(\epsilon_i), \sigma_{-1}(\mathbf{s}_{i_1}), \sigma_{-1}(\mathbf{s}_{i_2})$).

For the conjugated automorphism, since it involves a permutation of the polynomial's coefficients and a change of sign, it will suffice to check for each element that the permutation was done correctly as well as the sign change. In the case of the epsilon polynomial, one would only need to check the degree is at most $4$ by checking all other coefficients $\epsilon_{4} \dots \epsilon_{d-1}$ to be zero. For both cases, checking an element of a polynomial is equal to some element $b$ or to some element from another polynomial $c$ can be written as a dot product in this way: 

$$ct(\sigma_{-1}(X^{j})\mathbf{a} -b) = 0 \mod q'$$

$$ct(\sigma_{-1}(X^{j})\mathbf{a} -\sigma_{-1}(X^{k})\mathbf{c}) = 0 \mod q'$$







 

