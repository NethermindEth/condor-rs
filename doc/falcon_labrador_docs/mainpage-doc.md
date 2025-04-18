# Aggregating Falcon Signatures with LaBRADOR
This repository contains the implementation of *[Aggregating Falcon Signatures with LaBRADOR](https://eprint.iacr.org/2024/311.pdf)*.

The main goal of the paper—and of this implementation—is to aggregate Falcon signatures using a **non‑interactive** (via the Fiat–Shamir heuristic) version of LaBRADOR, a post‑quantum lattice‑based argument‑of‑knowledge scheme with short proofs.
The approach is to define a signature‑aggregation scheme (AS) for Falcon that relies on a LaBRADOR‑based succinct non‑interactive argument of knowledge (SNARK). 
In this SNARK we treat the signatures as **witnesses** and the messages and public keys as **statements**.  
The result is a non‑sequential aggregation scheme that requires zero interaction between signers, and, thanks to the succinct proofs, the scheme is well suited to bandwidth‑constrained settings such as blockchains.

These notes serve as a friendly introduction to the protocol and a prototype for the documentation. 
This implementation assumes Falcon signature scheme and interactive version of LaBRADOR are implemented.

## Overview 
The implemented changes to LaBRADOR consist of the following steps: 
- Changing the Modulus & Norm Checks
- Reformulating Constraints
- Working over Subring

## Assumption
We assume a user aims to aggregate $N$ signatures, with $i$ indexing each signature, message, and public key. 
For instance, $m\_i$ denotes the $i$-th message to be signed, $(\mathbf{s}\_{i\_1}, \mathbf{s}\_{i\_2})$ represents its signatures, and $pk\_i$ is the corresponding public key.


## Notation
Throughout this documentation let $q$ denote the modulus, and let $\mathbb{Z}_q$ be the ring of integers modulo $q$.  
Define the polynomial rings $\mathcal{R} = \mathbb{Z}[X]/(X^d + 1)$ and $\mathcal{R}_q = \mathbb{Z}_q[X]/(X^d + 1).$

We use the following conventions.
- **Non‑bold letters** (elements in $\mathbb{Z}_q$)  
  - $s \in \mathbb{Z}_q$: scalar  
  - $\vec{s} \in \mathbb{Z}_q^n$: vector of length $n$  
  - $A \in \mathbb{Z}_q^{m \times n}$: matrix with $m$ rows and $n$ columns  

- **Bold letters** (polynomial functions in $\mathcal{R}$ or $\mathcal{R}_q$)  
  - $\mathbf{s} \in \mathcal{R}$ or $\mathcal{R}_q$
  - $\vec{\mathbf{s}} \in (\mathcal{R}_q)^n$: vector of polynomials with length $n$ 
  - $\mathbf{A} \in (\mathcal{R}_q)^{m \times n}$: matrix of polynomials with $m$ rows and $n$ columns  

- $ct(\mathbf{f})$ denotes the constant term of the polynomial function $\mathbf{f} = a_0 + a_1X + \cdots + a_{d-1}X^{d-1},$ i.e. $ct(\mathbf{f}) = a_0$.
