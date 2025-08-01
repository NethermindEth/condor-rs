# Recusive Proof
To reduce the proof size, the prover **recursively** proves knowledge of the vectors instead of sending them outright.  

Concretely, the prover shows they know $\vec{\mathbf{z}},\\;\vec{\mathbf{t}},\\;\vec{\mathbf{g}},\\;\vec{\mathbf{h}}$ satisfying

$$
\vec{\mathbf{u}}\_{ 1}= \mathbf{B}\\,\vec{\mathbf{t}}
\\;\in\\; \mathcal{R}\_{ q}^{\kappa\_1},\qquad
\lVert\vec{\mathbf{t}}\rVert \\;\le\\; \gamma\_1.
\tag{1}
$$

$$
\mathbf{A}\\,\vec{\mathbf{z}}
\\;=\\;
\sum_{i= 1}^{ r}\mathbf{c}\_i\\,\vec{\mathbf{t}}\_i
\\;\in\\; \mathcal{R}\_{q}^{\kappa},\qquad
\lVert\vec{\mathbf{z}}\rVert \\;\le\\; \gamma.
\tag{2}
$$

$$
\langle\vec{\mathbf{z}},\vec{\mathbf{z}}\rangle=\sum_{i,j=1}^{r}\mathbf{g}\_{ij}\\,\mathbf{c}\_i\mathbf{c}\_j,\qquad
\sum_{i=1}^{r}
\langle\vec{\boldsymbol{\varphi}}\_i,\vec{\mathbf{z}}\rangle\\,\mathbf{c}\_i=\sum_{i,j=1}^{r}\mathbf{h}\_{ij}\\,\mathbf{c}\_i\mathbf{c}\_j,
$$

$$
\sum_{i,j= 1}^{ r}\mathbf{a}\_{ij}\\,\mathbf{g}\_{ij}
\\;+\\;
\sum_{i= 1}^{ r}\mathbf{h}\_{ii}
\\;-\\;
\mathbf{b}
\\;=\\;
\mathbf 0 .
\tag{3}
$$

$$
\vec{\mathbf{u}}\_{ 2}=\mathbf{C}\\,\vec{\mathbf{g}}
+
\mathbf{D}\\,\vec{\mathbf{h}}
\\;\in\\;
\mathcal{R}\_{q}^{\kappa\_{ 2}},\qquad
\sqrt{\lVert\vec{\mathbf{g}}\rVert^{ 2}
      +
      \lVert\vec{\mathbf{h}}\rVert^{ 2}}
\\;\le\\;
\gamma\_{ 2}.
\tag{4}
$$

We now show how to rewrite (**1**)–(**4**) as a new dot‑product instance that
LaBRADOR can handle.  The goal is to derive the public parameters and witness for the *next* recursion level.

Let  

$$
\vec{\mathbf{z}}=\vec{\mathbf{z}}^{(0)}
+
b\\,\vec{\mathbf{z}}^{( 1)},
\qquad\qquad
\vec{\mathbf{v}}=\vec{\mathbf{t}}\\;\Vert\\;\vec{\mathbf{g}}\\;\Vert\\;\vec{\mathbf{h}}
\\;\in\\;
\mathcal{R}\_{q}^{m}.
$$

The new witness is $\vec{\mathbf{z}}^{(0)} \\| \vec{\mathbf{z}}^{(1)} \\| \vec{\mathbf{t}} \\| \vec{\mathbf{g}} \\| \vec{\mathbf{h}}$. Following Section 5.3 of the paper we split this into $2v + u$ vectors (the new multiplicity) of rank $n'$.

The generic dot‑product constraint has the form

$$
\mathbf{g}^{(k)}\\!\bigl(\vec{\mathbf{s}}\_{\mathbf 1},\ldots,
                        \vec{\mathbf{s}}\_{\mathbf{r'}}\bigr)=\sum_{i,j= 1}^{r'}
      \mathbf{a}^{(k)}\_{ij}\\,
      \langle\\,\vec{\mathbf{s}}\_i,\vec{\mathbf{s}}\_j\rangle
+
\sum_{i=1}^{r'}
      \bigl\langle
      \vec{\boldsymbol{\varphi}}^{(k)}\_{i},
      \vec{\mathbf{s}}\_i
      \bigr\rangle-\mathbf{b}^{(k)}=\mathbf 0 ,
\tag{6}
$$

for $k= 1,\ldots,\kappa
+
\kappa_{1}
+
\kappa_{2} + 3=
K'$. 
Assume each vector $\vec{\boldsymbol{\varphi}}_i^{(k)}$ has length $m_1$, and there are $r'$ such vectors. Concatenate all their coordinates into one long vector and denote its entries by $\boldsymbol{\varphi}_i^{(k)}$ for $i \in [1, m_1 r']$. With this flattening, Equation (6) becomes a single inner-product of the two concatenated vectors, instead of a sum over $r'$ smaller inner-products. 


To satisfy (1)-(4), we need  
1.  3 constraints for (**3**)
2. $\kappa$ for the vector equation (**2**), and similarly $\kappa_1$ and $\kappa_2$ for (**1**) and (**4**).


## Equations

### Equation (1)
*(to be filled)*

### Equation (2)
We want to cast

$$
\mathbf{A}\\,\vec{\mathbf{z}}=
\sum_{i=1}^{r}\mathbf{c}\_i\\,\vec{\mathbf{t}}\_i,
\qquad
\mathbf{c}\_i\in\mathcal{R}\_{q},\\;
\mathbf{A}\vec{\mathbf{z}}\in\mathcal{R}\_{q}^{\kappa},
$$

into the dot‑product template. 

Each inner commitment is of the form:

$$
\vec{\mathbf{t}}\_i=
\vec{\mathbf{t}}^{(0)}\_i
+
b\_{1}\\,\vec{\mathbf{t}}^{(1)}\_i
+
\ldots
+
b\_{1}^{t\_{1}- 1}\\,
\vec{\mathbf{t}}^{(t\_{1}-1)}\_i .
$$

Hence  

$$
\mathbf{A}\\,\vec{\mathbf{z}}^{(0)}+
b\\,\mathbf{A}\\,\vec{\mathbf{z}}^{(1)}-
\sum_{i=1}^{r}
\Bigl(
\mathbf{c}\_i\vec{\mathbf{t}}^{(0)}\_i+
b\_{1}\mathbf{c}\_i\vec{\mathbf{t}}^{(1)}\_i+
\ldots+
b\_{1}^{t\_{1}-1}
\mathbf{c}\_i\vec{\mathbf{t}}^{(t\_{1}-1)}\_i
\Bigr)=
\mathbf{\vec{0}}.
$$

Because the left‑hand side is a length‑$`\kappa`$ vector, we use $` \kappa`$ separate constraints—one per coordinate $`k \in [1 , \kappa]`$. For each coordinate $`k`$, we have:

$$
\langle\vec{\mathbf{A}}\_k,\vec{\mathbf{z}}^{(0)}\rangle+
\langle b\\,\vec{\mathbf{A}}\_k,\vec{\mathbf{z}}^{(1)}\rangle-
\sum_{i=1}^{r}
\Bigl(
\mathbf{c}\_i{\mathbf{t}}^{(0)}\_{ik}+
b\_{1}\mathbf{c}\_i{\mathbf{t}}^{(1)}\_{ik}+
\ldots+
b\_{1}^{t\_{1}-1}
\mathbf{c}\_i{\mathbf{t}}^{(t\_{1}-1)}\_{ik}
\Bigr)=
\mathbf 0.
$$ 

where $`\vec{\mathbf{A}}_k`$ is the $`k`$-th row of matrix $`\mathbf{A}`$, and $`\mathbf{t}^{(j)}_{ik}`$ is the $`k`$-th row of vector $`\vec{\mathbf{t}}^{(j)}_{i}`$. If we concatinate all vectors $`\vec{\mathbf{t}}^{(j)}_{i}`$ for all $`i`$ and $`j`$ to one vector $`\mathbf{\vec{t}}`$, element $`\mathbf{t}^{(j)}_{ik}`$ is the 
$`(i-1) t_{ 1}\kappa+j\kappa + k`$-th element of the vector $`\mathbf{\vec{t}}`$.

Now, we can cast the above equation to dot product relation as follow:

$$
\mathbf{a}\_{ij}^{(k)} = \mathbf 0, \qquad
\boldsymbol{\varphi}\_{i}^{(k)} = \mathbf{A}\_{ki}   \qquad i\in[1,n], \qquad \boldsymbol{\varphi}\_{i}^{(k)} = \mathbf{A}\_{k(i-n)}  \qquad i\in[ n\\!+\\! 1, 2n]
$$

$$
\boldsymbol{\varphi}_{\\,2n+(i-1) t\_{ 1}\kappa+j\kappa + k}^{(k)}=
\mathbf{c}\_i\\,b\_{1}^{\\,j},
\quad
i\in[ 1, r],\\;
j\in[ 0,t\_{ 1}\\!-\\! 1],
$$

and all remaining public coefficients are set to **zero**.

### Equation (3)

### Equation (4)

## Norm Check

