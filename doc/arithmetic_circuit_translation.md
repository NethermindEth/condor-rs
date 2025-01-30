## Step 1 - Arithmetic Circuit Translation

In the world of computational complexity, one of the foundational results is the Cook-Levin theorem, which established the concept of NP-completeness. This theorem demonstrated that a broad range of problems in the NP class could be reduced to one central problem: **Boolean satisfiability (SAT)**. This means that it allows us to express common NP problems instances as a set of Boolean expressions that must be satisfied. These problems, typically modeled using Boolean logic, can also be rewritten as arithmetic circuits. These circuits, which use basic operations like addition and multiplication, can represent NP-complete problems in a new way.

The main idea regarding this proof systems involving circuits, its to write a statement one wants to prove (normally that some calculation has been done) by rewriting that calculation as an arithmetic circuit or an equivalent formulation that will allow for automatic proof protocols to convince a verifier of such statement.

### Definition of an Arithmetic Circuit

The main idea behind arithmetic circuits is that, given an instance of a hard problem, we can rewrite it as a **Directed Acyclic Graph (DAG)**, which can be imagined as a tree for simplicity. Each node in the graph represents either a sum or a multiplication, and the leaves of the tree represent the variables and constants in question. Every NP problem can be rewritten in this form, much like how every piece of code can be converted into a circuit of binary operations.

For example, consider the equation $x + y$, where both $x$ and $y$ are variables. This can be represented as a tree with three nodes: a central node representing addition and two leaf nodes, each connected to the central node, representing the variables $x$ and $y$.

This kind of arithmetic circuit $C$ allows us to pose problems by reducing them to **arithmetic constraints**. A solution to such a problem would be an $s$ such that $C(s) = 0$. It’s important to note that the problem we define here involves **multivariate polynomials**, since all the operations ultimately form polynomials based on the variables at the leaves of the DAG.

### R1CS and Dot Product Constraints

Arithmetic circuits can be represented in a standard way called a **Rank-1 Constraint System (R1CS)**, where given three matrices $A, B, C$, an element $s$ all modulo $N$ is a solution if:

$$
(A s) \circ (B s) = C s \text{ mod } N
$$

where $\circ$ denotes the **Hadamard product** (element-wise multiplication), and the rest follows standard matrix-vector multiplication.  

A simple way to observe that an arithmetic circuit can be represented in this form is to understand that an R1CS system consists solely of **multiplications and additions** that generate constraints. By writing:

$$
(As \circ Bs) - Cs = 0 \text{ mod } N
$$

and rewriting the right-hand side as a **Directed Acyclic Graph (DAG)**, as explained earlier, the solutions $s$ correspond to the values assigned to the leaves, ensuring that the circuit evaluates to zero.

### Dot Product Constraints in LaBRADOR

**LaBRADOR** was specifically designed to work with **dot product constraints**, which are constraints of the form:

$$
f(s)=0 \text{ or } ct(f(s)) = 0
$$

where $ct()$ extracts the constant coefficient of the polynomial. The function $f(s)$ has the general form:

$$
f(s) = \sum_{1 \leq i,j \leq r} a_{ij} \langle s_i, s_j \rangle + \sum_{i=1}^{r} \langle \varphi_i, s_i \rangle + b
$$

where both the variables $s$ and the constants $a_{ij}, \varphi_i, b$ are **polynomials (or vectors of polynomials) over a cyclotomic ring**:

$$
Z_{q}\[x\] / (x^d + 1)
$$

This structure is particularly useful in **lattice-based cryptography** and **zero-knowledge proof systems**.


### Binary R1CS to Dot Product Constrains

### R1CS modulo $2^{d}+1$ to Dot Product Constrains
