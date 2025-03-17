# Aggregating Falcon Signatures with LaBRADOR

This is the code implementation of "Aggregating Falcon Signatures with LaBRADOR." A non-interactive version of LaBRADOR, utilizing the Fiat-Shamir heuristic, that allows for a significant reduction in proof sizes within a signature aggregation scheme, compared to a basic concatenation procedure, while still relying on the security of standard lattice problems.

These notes serve as a friendly introduction to the protocol and a prototype for the documentation. They are based on the assumption that one has already implemented Falcon and the original interactive version of LaBRADOR.

The main idea is to define a signature aggregation scheme (AS) for Falcon based on the use of a succinct non-interactive argument of knowledge (SNARK), where we can set the signatures as witnesses and the messages and public keys as statements. This would allow for a non-sequential signature aggregation scheme, allowing zero interaction between signers. Additionally, the succinctness of the SNARK is perfect for bandwidth bottleneck situations, such as in blockchain. This code is an adaptation of LaBRADOR to function with Falcon for AS.

## Overview 

The implemented changes consist of the following steps: 
- Changing the Modulus & Norm Checks
- Reformulating Constraints
- Working over Subring




## Notation
We will use an upper bar $\bar{s}$ for vectors, lowercase $s$ for scalars, uppercase $S$ for matrices, and boldface letters for elements $\mathbf{s} \in \mathbb{Z}_q\[x\] / (x^d + 1)$, unless explicitly noted.
