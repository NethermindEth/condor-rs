## Step 5 - Amortization
Just like in the projection step, were we used a statistical argument to verify that that the norm of the witness vector is small, in the Amortization step both the correctness of the second step in aggregation and the verification of the commitments are proven probabilistically. 

In this step the verifier sends a colection of random polynomials $\mathbf{c_{1}},\dots,\mathbf{c_{r}}$ as challenges, therefore defining a random linear combination of openings $\mathbf{\bar{z}}=\mathbf{c_{1}\bar{s_{1}}}+\dots+\mathbf{c_{r}\bar{s_{r}}}$.
In order to verify that the Ajtai commitments where done correctly, instead of sending the wintness vectors the verifyer can use $\mathbf{\bar{z}}$ to gain information on the witness vectors, more specifically, the verifier can check both that the norm of $\mathbf{\bar{z}}$ is suficiently small and also that $A\mathbf{\bar{z}} = \sum_{i=1}\^{r}\mathbf{c_{i}}A\mathbf{\bar{s_{i}}} =\sum_{i=1}\^{r}\mathbf{c_{i}}\mathbf{\bar{t_{i}}}$ given that it is in posession of $A$, $\mathbf{\bar{z}}$, $t$ and the challenges $\mathbf{c_{i}}$.


In the last aggregation step, the verifier sends uniformly distributed challenges $\bar{\alpha} \in (Z_{q})\^{L}$ and $\bar{\beta} \in (Z_{q})\^{\lceil 128/ log(q) \rceil}$,such that the same scheme could be applied again, but now for all the dot product constraints (including those sent in the first aggregation step).

So, we will aggregate all functions into one linear combination of functions of the form:
$$\sum_{k=1}\^{K}\bar{\alpha_{k}}f\^{\(k\)}(\mathbf{\bar{s_{1}}},\dots,\mathbf{\bar{s_{r}}})+\sum_{k=1}\^{\lceil 128/ log(q) \rceil}\bar{\beta_{k}}f^{''\(k\)}(\mathbf{\bar{s_{1}}},\dots,\mathbf{\bar{s_{r}}})$$

When we decompose these functions into their respective parts, we will only be interested in sending  $\mathbf{h_{ij}} = \frac{1}{2}(<\bar{\varphi_{i}},\mathbf{\bar{s_{j}}}>+<\mathbf{\varphi_{j}},\mathbf{\bar{s_{i}}}>)$ since the information regarding $<\mathbf{\bar{s_{i}}},\mathbf{\bar{s_{j}}}>$ was already sent with the polynomial $\mathbf{\bar{g}}$ on the Ajtai commitment step.

The $\bar{\varphi_{i}}$ that were mentioned in $\mathbf{h_{ij}}$ are defined as:
$$\sum_{k=1}\^{K}\bar{\alpha_{k}}\bar{\varphi_{i}}\^{\(k\)}+\sum_{k=1}\^{\lceil 128/ log(q) \rceil}\bar{\beta_{k}}\bar{\varphi_{i}}\^{''\(k\)}$$
Just like in the case of the garbage polynomail $\mathbf{\bar{g}}$ the $\mathbf{h_{ij}}$ also form a matrix of polynomials that can be decomposed based on some base $b_{1}$ and sent in the same way as $\mathbf{\bar{g}}$ by defining: 
$$\mathbf{\bar{u_{2}}}=D*\mathbf{\bar{h}}$$

