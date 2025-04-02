# Reformulating Constraints
In order to achieve an improvement in computational complexity, witnesses will be rewritten via a padding scheme. This new format, although increasing both in size and in the number of dot product constraints that are required, will significantly reduce the need for garbage polynomials, thereby providing a general increase in speed.

## Padding Scheme
In this padding scheme, witness elements $w_{1} \dots w_{N}$ $\in \mathcal{R}_{q\^{'}}$ along with other elements $w_{1}^{'} \dots w_{N}^{'} \in \mathcal{R}_{q\^{'}}$ are such that $\lVert w_{i} \rVert_{2}^{2} = ct(w_{i}' w_{i})$

(An example of $w_{i}'$ being the conjugate automorphism of $w_{i}$) can be rewritten as vectors of the form $\vec{u}_{1}, \dots, \vec{u}_{\lceil N/\rho \rceil} \in \mathcal{R}_{q'}^{N}$ and $\vec{u}_{1}', \dots, \vec{u}_{\lceil N/\rho \rceil}' \in \mathcal{R}_{q'}^{N}$ respectively, where $\rho = \lfloor \sqrt{N} \rfloor$.

After following the padding equations described in the paper, we can think of each vector $\vec{u}_{i}$ as a vector filled with zeros except for some locations where the original witness polynomials will be located: $[0, 0, 0, \dots, w_{2\rho +1}, \dots, w_{3\rho}, \dots, 0, 0]$. We can access each old witness polynomial using the index functions of the form:
$$ index(i) = \lceil i/\rho \rceil $$
$$ index'(i) = ((i-1) \mod \rho) + 1 $$

for $u_{i}$ and $u_{i}'$ respectively. From this, it follows that we can write $\lVert w_{i} \rVert_{2}^{2} = ct(w_{i}'w_{i}) = ct(\langle \vec{u}_{index'(i)}', \vec{u}_{index(i)} \rangle)$.

Beyond adding dot product restrictions to check each zero element of the new witness vectors, we will also rewrite all the restrictions using the new witnesses as well as the index function. In order to rewrite the constraints, there will be extensive use of $\delta_{i} \in \mathcal{R}_{q'}^{N}$, the vector with the $i$-th entry as $1$ (the identity of the ring) and all other elements as $0$. 
As an example, we can rewrite Falcon's restrictions as:
$$ \langle \delta_{i}, \vec{u}_{index(i)} \rangle + \langle h_{i} \delta_{i}, \vec{u}_{index(i)} \rangle + \langle q \delta_{i}, v \rangle - t_{i} = 0 $$

