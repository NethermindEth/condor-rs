## Step 2 - Ajtai Commitment
In the first step of the protocol, the prover commits to the small-norm solution vectors  $s_{1},...,s_{r}$ by computing **Ajtai commitments**. 
The idea behind a commitment is a cryptographic equivalent to keeping some knowledge sealed in an **'envelope'**, to be revealed later.
Specifically, the prover computes the vector $t_{i} = A * s_{i}$, where $s_{i}$ is a vector of $n$ polynomials. 
After multiplication, this results in vectors $t_{i}$ of $k$ polynomials. 
However, sending the vectors $t_{i}$ directly could be costly due to its potentially large size. 
Therefore, the prover commits to the $t_{i}$ using another Ajtai commitment, referred to as the **outer commitment**.

### Decomposition for Efficiency

The components of the vector $t_(i)$ have coefficients that are arbitrary modulo $q$. 
This may lead to inefficiencies when transferring large numbers. To optimize the transfer and reduce the size of the values being committed, 
the coefficients need to be decomposed into smaller parts.

The decomposition is done **element-wise**. Specifically, each polynomial in $t_{i}$ is decomposed into $t_{1} \geq 2$ parts with respect to a small base $b_{1}$, 
such that:

$$t_{i} = t_{i}\^{\(0\)} + t_{i}\^{1} * b_{1} + t_{i}\^{\(2\)} * b_{1}^{2} + ... + t_{i}^{\(t_{1}-1\)} * b_{1}^{t_{1}-1}$$

In this decomposition, centered representatives are used, which ensures that each element of the vector lies within the range $\[\frac{-b_{1}}{2}, \frac{b_{1}}{2}\]$.
Once the decomposition is complete for each $t_{i}$, we concatenate all the decomposition coefficients $t_{i}\^{k}$ for each $i$ and $k$ to form a new vector $t$. 
The second Ajtai commitment can then be written as:
$u_{1} = B * t$

### Commitment to the Garbage Polynomial

The protocol also includes a polynomial referred to as the garbage polynomial, which does not depend on any challenge during the interaction and is used in the amortization part. To simplify the proof of security, this polynomial is also committed to at the beginning of the protocol. To commit to the garbage polynomial, the prover can simply include it in the commitment to $t$. 

The garbage polynomial commitment is handled similarly. Given that $g_{ij} = <s_{i}, s_{j}>$ is the inner product of the solution vectors, 
the prover constructs a symmetric matrix of polynomials. Each element of this matrix $g_{ij}$ is decomposed based on some base $b_{2}$ into $t_{2} \geq 2$ parts and then each decomposition coefficient is concateated into $g$ :

$$g_{ij} = g_{ij}\^{\(0\)} + g_{ij}^{\(1\)} * b_{2} + g_{ij}\^{\(2\)} * b_{2}\^{2} + ... + g_{ij}\^{\(t_{2}-1\)} * b_{2}\^{t_2-1}$$

Finally, the full commitment to both the vector $t$ and the garbage polynomial $g$ is expressed as:
$$u_{1} = B * t + C * g$$

### Norm Constraints & Common Reference String
- $A$, $B$, and $C$ are matrices that serve as common references between the verifier and the prover
- ||$t$|| $\leq γ_{1}$
- ||$g$|| $\leq γ_{2}$

(More details on choosing the bases and $γ$ will be added soon)

