## Step 5 - Amortization
Just like in the projection step, were we used a statistical argument to verify that that the norm of the witness vector is small, in the Amortization step both the correctness of the second step in aggregation and the verification of the commitments are proven probabilistically. 

In this step the verifier sends a colection of random polynomials $c_{1},..,c_{r}$ as challenges, therefor defining a random linear combination of openings $z=c_{1}t_{1}+...+c_{r}t_{r}$.
In order to verify that the Ajtai commitments where done correctly, instead of sending the wintness vectors the verifyer can use $z$ to gain information on the witness vectors, more specifically, the verifier can check both that the norm of $z$ is suficiently small and also that $Az = \sum_{i=1}\^{r}c_{i}t_{i}$


In the last aggregation step, the verifier sends uniformly distributed challenges $\alpha \in (Z_{q})\^{L}$ and $\beta \in (Z_{q})\^{\lceil 128/ log(q) \rceil} $,such that the same scheme could be applied again, but now for all the dot product constraints (including those sent in the first aggregation step).

So, the functions we want to send are of the form:
$$\sum_{k=1}\^{K}\alpha_{k}f\^{\(k\)}(s_{1},...,s_{r})+\sum_{k=1}\^{\lceil 128/ log(q) \rceil}\beta_{k}f^{''\(k\)}(s_{1},...,s_{r})$$
When we decompose these functions into their respective parts, we will only be interested in sending  $h_{ij} = \frac{1}{2}(<\varphi_{i},s_{j}>+<\varphi_{j},s_{i}>)$ since the information regarding $<s_{i},s_{j}>$ was already sent with the polynomial $g$ on the Ajtai commitment step.

The $\varphi_{i}$ that where mentioned in $h_{ij}$ are defined as:
$$\sum_{k=1}\^{K}\alpha_{k}\varphi_{i}\^{\(k\)}+\sum_{k=1}\^{\lceil 128/ log(q) \rceil}\beta_{k}\varphi_{i}\^{''\(k\)}$$
Just like in the case of the garbage polynomail $g$ the $h_{ij}$ also form a matrix of polynomials that can be decomposed based on some base $b_{1}$ and sent in the same way as $g$ by defining: 
$$u_{2} = D*h$$

