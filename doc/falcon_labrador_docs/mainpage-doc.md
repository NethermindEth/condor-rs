# Aggregating Falcon Signatures with LaBRADOR

This is the code implementation of "Aggregating Falcon Signatures with LaBRADOR." A non-interactive version of LaBRADOR, utilizing the Fiat-Shamir heuristic, that allows for a significant reduction in proof sizes within a signature aggregation scheme, compared to a basic concatenation procedure, while still relying on the security of standard lattice problems.

These notes serve as a friendly introduction to the protocol and a prototype for the documentation. They are based on the assumption that one has already implemented Falcon and the original interactive version of LaBRADOR.

The main idea is to define a protocol between a prover $P$ and a verifier $V$, where, given an initial arithmetic circuit, they can perform an interaction $<P, V>$. This interaction allows the prover to demonstrate knowledge of solutions to the circuit, and the verifier will output a binary result indicating whether the proof was accepted.  
Protocols of this type already exist, but LaBRADOR has the added advantage of achieving sublinear proof sizes while relying on (plausible) post-quantum secure cryptographic lattice problems.

## Protocol Overview

The protocol, from the circuit to the final proof, can be described through the following sequence of steps:  
1. Translate the rank-1 constraint system (R1CS) into a dot product constraint system and solve it.  
2. Apply an Ajtai commitment to the solutions of the dot product constraints.  
3. Prove the norm requirements on the solution using projections (Modular Johnson-Lindenstrauss Lemma).  
4. Make the proof more compact by aggregating multiple constraints (Aggregation).  
5. Only prove linear combinations of the commitments (Amortization).  Instead of proving all commitments directly, prove linear combinations of them. This recursive amortization process reduces the witness size in each iteration.


## Notation
We will use an upper bar $\bar{s}$ for vectors, lowercase $s$ for scalars, uppercase $S$ for matrices, and boldface letters for elements $\mathbf{s} \in \mathbb{Z}_q\[x\] / (x^d + 1)$, unless explicitly noted.
