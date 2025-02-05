## Step 2 - Ajtai Commitment

Many well-known computational problems, such as factoring, are considered worst-case problems, whereas cryptographic applications typically require problems that are hard on average. This is because, for hard average-case problems, any random instance is likely to be hard with a positive probability. In contrast, there is no known method for generating a hard instance of a worst-case problem. By leveraging lattice-based problems, it's possible to construct average-case problems that are as difficult as certain well-known worst-case problems.

In 1996, Miklós Ajtai proposed the idea of how to generate one of these average hard instance problems. An Ajtai Commitment is a polynomial commitment scheme based on these techniques to guarantee the average hardness of this fundamental component in many cryptographic protocols.

### What is an Ajtai Commitment?

An Ajtai commitment scheme allows a committer to publish a value, called the commitment, which binds them to a message without revealing it. Later, they may open the commitment and reveal the committed message to a verifier, who can check that the message is consistent with the commitment. In our case, we are more interested in the compactness of the commitment in relation to the message, allowing the committer to send smaller messages, rather than the non-revealing aspect. The idea behind a commitment is a cryptographic equivalent to keeping some knowledge sealed in an **"envelope"**, to be revealed later.

The fundamental problem upon which it is based is called the **Modular-SIS** problem. The Short Integer Solution (SIS) problem is based on the idea that, given a sufficiently large random matrix $A$ and a vector $\bar{x}$, solving the linear system of equations $A\bar{x}=\bar{b}$, where the solution $\textbf{x}$ is required to have a small norm $\lVert \textbf{x} \rVert < \beta$, is as hard as the **short basis problem** (a lattice problem that is believed to be hard on average). The **Modular** part is simply a generalization of SIS to work on a module structure.

In order to perform an Ajtai Commitment Scheme, the committer would first need to be in possession of a **small-norm** vector $\bar{x}$ that represents the message they want to commit to. Given a public matrix $A$, one can calculate $A \bar{x} = \bar{b}$ and send $\bar{b}$, which is a more compact vector, making it lighter to send over restricted-size networks. Now, the verifier is in possession of both $A$ and $\bar{b}$ (the "envelope"). Once the solution vector $\bar{x}$ is revealed at the end by the prover, the verifier can easily check that the prover had the solution all along.

The reason why $\bar{b}$ is binding is due to the SIS problem. It is a hard problem to find a solution to the constraints that is also of small norm, so the prover must have had the solution all along.


### Commitment
The prover in the protocol has solved the dot product constraints and is in possession of small-norm solution vectors (which we will call "witness vectors") as well as public matrices $A,B,C,D$.

In the first step of the protocol, the prover commits to the small-norm solution vectors  $\mathbf{\bar{s_{1}}},\dots,\mathbf{\bar{s_{r}}}$ by computing **Ajtai commitments**. 
Specifically, the prover computes the vector $\mathbf{\bar{t_{i}}} = A * \mathbf{\bar{s_{i}}}$, where $\mathbf{\bar{s_{i}}}$ is a **vector** of $n$ polynomials. 
After multiplication, this results in vectors $\mathbf{\bar{t_{i}}}$ of $k$ polynomials. 
However, sending the vectors $\mathbf{\bar{t_{i}}}$ directly could be costly due to its potentially large size. 
Therefore, the prover commits to the $\mathbf{\bar{t_{i}}}$ using another Ajtai commitment, referred to as the **outer commitment**.

### Decomposition for Efficiency

The components of the vector $\mathbf{\bar{t_{i}}}$ have coefficients that are arbitrary modulo $q$. 
This may lead to inefficiencies when transferring large numbers. To optimize the transfer and reduce the size of the values being committed, 
the coefficients need to be decomposed into smaller parts.

The decomposition is done **element-wise**. Specifically, each polynomial in $\mathbf{\bar{t_{i}}}$ is decomposed into $t_{1} \geq 2$ parts with respect to a small base $b_{1}$, 
such that:

$$\mathbf{\bar{t_{i}}} = \mathbf{\bar{t_{i}}}\^{\(0\)} + \mathbf{\bar{t_{i}}}\^{\(1\)} * b_{1} + \mathbf{\bar{t_{i}}}\^{\(2\)} * b_{1}^{2} + \dots + \mathbf{\bar{t_{i}}}^{\(t_{1}-1\)} * b_{1}^{t_{1}-1}$$

In this decomposition, centered representatives modulo $b_{1}$ are used, which means $\lVert \mathbf{\bar{t_{i}}}\^{\(k\)} \rVert_{\infty}$
Once the decomposition is complete for each $\mathbf{\bar{t_{i}}}$, we concatenate all the decomposition coefficients $\mathbf{\bar{t_{i}}}\^{k}$ for each $i$ and $k$ to form a new vector $\mathbf{\bar{t}}$. 
The second Ajtai commitment can then be written as:
$\mathbf{\bar{u_{1}}}=B*\mathbf{\bar{t}}$



### Commitment to the Garbage Polynomial

The protocol also includes a polynomial referred to as the garbage polynomial, which does not depend on any challenge during the interaction and is used in the amortization part. To simplify the proof of security, this polynomial is also committed to at the beginning of the protocol. To commit to the garbage polynomial, the prover can simply include it in the commitment to $\mathbf{\bar{t}}$. 

The garbage polynomial commitment is handled similarly. Given that $\mathbf{g_{ij}} = <\mathbf{\bar{s_{i}}}, \mathbf{\bar{s_{j}}}>$ is the inner product of the solution vectors, 
the prover constructs a symmetric matrix of polynomials. Each element of this matrix $\mathbf{g_{ij}}$ is decomposed based on some base $b_{2}$ into $t_{2} \geq 2$ parts and then each decomposition coefficient is concateated into $\mathbf{g}$ :

$$\mathbf{g_{ij}} = \mathbf{g_{ij}}\^{\(0\)} + \mathbf{g_{ij}}\^{\(1\)} * b_{2} + \mathbf{g_{ij}}\^{\(2\)} * b_{2}\^{2} + \dots + \mathbf{g_{ij}}\^{\(t_{2}-1\)} * b_{2}\^{t_2-1}$$

Finally, the full commitment to both the vector $\mathbf{\bar{t}}$ and the garbage polynomial $\mathbf{\bar{g}}$ is expressed as:
$$\mathbf{\bar{u_{1}}} = B * \mathbf{\bar{t}} + C * \mathbf{\bar{g}}$$


