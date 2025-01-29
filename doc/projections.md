## Step 3 - Projections (Johnson-Lindenstrauss)
The main idea behind the Modular Johnson-Lindenstrauss Lemma is that if you want to prove that a long vector has a small norm, instead of having to send the entire vector to verify this condition, the verifier can send random linear projections $\Pi_{i} $ as challenges to the prover and then the prover just needs to return the aplied proyection to the vector $p = \sum_{i=1}\^{r}\Pi s_{i}$ which is much more compact. The lemma guarantees that the projections almost preserve the L2 norm, so by performing a sufficient number of challenges, one can gain information about the actual norm size of the original vector.

After the projection vectors are sent to the prover, so that it can check that they have small norms, a proof is also required to verify that the projections were actually performed correctly. To do this, the prover will transform the projections into dot product constraints. These constraints will later be sent as part of the aggregation step. The way this works is by rewriting this matrix multiplication $\Pi_{i}s_{i}$ as a sum of dot products $\sum_{i=1}^{r} \langle \pi_{i}^{(j)}, s_{i} \rangle$, Where this is the inner product between the coefficients of the polynomials and the j-th row of the projection matrix. Finally, we have the projection vector $p = (p_{1}, ..., p_{256})$, where each coordinate can be rewritten as $p_{j} = \sum_{i=1}^{r} \langle \pi_{i}^{(j)}, s_{i} \rangle$, and this allows us to define the dot product constraints:

$$ \sum_{i=1}^{r} \langle \pi_{i}^{(j)}, X \rangle - p_j = 0 \quad \forall j \in \{1, \dots, 256\} $$

So that this dot product constrain becomes zero at $s_{1}, \dots, s_{r}$.

