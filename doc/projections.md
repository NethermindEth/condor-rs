## Step 3 - Projections (Johnson-Lindenstrauss)
The main idea behind the Modular Johnson-Lindenstrauss Lemma is that if you want to prove that a long vector has a small norm, instead of having to send the entire vector to verify this condition, the verifier can send random linear projections $\Pi_{i}$ as challenges to the prover and then the prover just needs to return the aplied proyection to the vector $\bar{p} = \sum_{i=1}\^{r}\Pi \bar{s_{i}}$ which is much more compact. The lemma guarantees that the projections almost preserve the L2 norm, so by performing a sufficient number of challenges, one can gain information about the actual norm size of the original vector.

After the projection vectors are sent to the prover, so that it can check that they have small norms, a proof is also required to verify that the projections were actually performed correctly. To do this, the prover will transform the projections into dot product constraints. These constraints will later be sent as part of the aggregation step. The way this works is by rewriting this matrix multiplication $\Pi_{i}\bar{s_{i}}$ as a sum of dot products $\sum_{i=1}^{r} \langle \bar{\pi_{i}}^{(j)}, \bar{s_{i}} \rangle$, Where this is the inner product between the coefficients of the polynomials and the j-th row of the projection matrix. Finally, we have the projection vector $\bar{p} = (p_{1}, \dots, p_{256})$, where each coordinate can be rewritten as $p_{j} = \sum_{i=1}^{r} \langle \bar{\pi_{i}}^{(j)}, \bar{s_{i}} \rangle$, and this allows us to define the dot product constraints:

$$\sum_{i=1}^{r} \langle \sigma_{-1}\(\mathbf{\pi_{i}^{(j)}}\), \mathbf{X} \rangle - p_j  \quad \forall j \in \{1, \dots, 256\}$$

Let's explain why this makes sense. The expression $p_j = \sum_{i=1}^{r} \langle \bar{\pi_{i}}^{(j)}, \bar{s_{i}} \rangle$ is scalar-valued, so if we were to write $\sum_{i=1}^{r} \langle \bar{\pi_{i}}^{(j)}, \bar{s_{i}} \rangle - p_j$, the result would always be zero. To avoid this triviality and prove correctness, we need to rewrite it as a dot product constraint.

To accomplish this, we treat the values $\bar{\pi_{i}}^{(j)}$ as the coefficients of a polynomial. This is possible due to the one-to-one correspondence between a degree $d$ polynomial and a vector of length $d + 1$. Now, we can use the same polynomial in a dot product constraint. Specifically, we will utilize the automorphism $\sigma_{-1}$, defined as $\sigma_{-1}: X \to X^{-1}$. Applying this automorphism in our inner product ensures that the constant coefficient is $\sum_{i=1}^{r} \langle \bar{\pi_{i}}^{(j)}, \bar{s_{i}} \rangle$. Thus, we obtain the following equality:

$$
ct\left( \sum_{i=1}^{r} \langle \sigma_{-1}(\mathbf{\pi_{i}^{(j)}}), \mathbf{X} \rangle \right) - p_j = 0 \quad \forall j \in \{1, \dots, 256\}
$$

So the previous dot product belongs to the family of dot products with a zero constant coefficient under the variables $\mathbf{\bar{s_1}}, \dots, \mathbf{\bar{s_r}}$ and can be aggregated with other functions of the same form.


