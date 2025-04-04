# Reformulating Constraints
In order to achieve an improvement in computational complexity, witnesses will be rewritten via a padding scheme. This new format, although increasing both in size and in the number of dot product constraints that are required, will significantly reduce the need for garbage polynomials, thereby providing a general increase in speed.

## Padding Scheme
In this padding scheme, witness elements $w\_{1} \dots w\_{N}$ $\in \mathcal{R}\_{q\^{'}}$ along with other elements $w\_{1}\^{'} \dots w\_{N}\^{'} \in \mathcal{R}\_{q\^{'}}$ are such that $\lVert w\_{i} \rVert\_{2}\^{2} = ct(w\_{i}' w\_{i})$ (An example of $w\_{i}'$ being the conjugate automorphism of $w\_{i}$) can be rewritten as vectors of the form $\vec{u}\_{1}, \dots, \vec{u}\_{\lceil N \rho \rceil} \in \mathcal{R}\_{q'}\^{N}$ and $\vec{u}\_{1}', \dots, \vec{u}\_{\lceil N \rho \rceil}' \in \mathcal{R}\_{q'}\^{N}$ respectively, where $\rho = \lfloor \sqrt{N} \rfloor$.

After following the padding equations described in the paper, we can think of each vector $\vec{u}\_{i}$ as a vector filled with zeros except for some locations where the original witness polynomials will be located: $\[0, 0, 0, \dots, w_{2\rho +1}, \dots, w_{3\rho}, \dots, 0, 0\]$. We can access each old witness polynomial using the index functions of the form:
$$index(i) = \lceil i \rho \rceil$$
$$index'(i) = ((i-1) \mod \rho) + 1$$

for $u\_{i}$ and $u\_{i}'$ respectively. From this, it follows that we can write $\lVert w\_{i} \rVert\_{2}\^{2} = ct(w\_{i}'w\_{i}) = ct(\langle \vec{u}\_{index'(i)}', \vec{u}\_{index(i)} \rangle)$.

Beyond adding dot product restrictions to check each zero element of the new witness vectors, we will also rewrite all the restrictions using the new witnesses as well as the index function. In order to rewrite the constraints, there will be extensive use of $\delta\_{i} \in \mathcal{R}\_{q'}\^{N}$, the vector with the $i$-th entry as $1$ (the identity of the ring) and all other elements as $0$. 
As an example, we can rewrite Falcon's restrictions as:
$$\langle \delta\_{i}, \vec{u}\_{index(i)} \rangle + \langle h\_{i} \delta\_{i}, \vec{u}\_{index(i)} \rangle + \langle q \delta\_{i}, v \rangle - t\_{i} = 0$$

