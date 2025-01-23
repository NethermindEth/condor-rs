# LaBRADOR: Lattice-Based Recursively Amortized Demonstration Of R1CS

This is the code implementation of LaBRADOR, a Lattice based proof system for R1CS that can achive proof sizes of logarithmic order.
This notes will serve as a friendly introduction to the protocol and as a documentation prototype.
The main idea is to define a protocol between a prover P and a verifier V, that given some initial arithmetic circuit they can perform an interaction <P,V> that would allow the prover to prove knowledge of solutions to that circuit to the verifier, such that the verifier will output a binary output representing whether the proof was accepted or not. Protocols of this type already exist yet this one has the added advantage of achiving sublinear proof sizes while relying on (Provable) post quantum secure cryptographic lattice problems. 

(This documentation is under construction)

## Protocol
We can think of the protocol, from circuit all the way to final proof as this sequence of steps: 
    1 - Translate rank-1 constrain system to a dot product constrain system and solve it
    2 - Ajtai commitment to the solutions of the dot product constrains (witness vectors)
    3 - Prove norm requirments on solution using proyections (Modular Johnson-Lindenstrauss Lemma) 
    4 - Make the proof more compact by aggregating multiple constrains (Aggregetion)
    5 - Only prove linear combinations of the commitments (Amortization)
    
## Step 1 - Artithmetic circuit translation 

## Step 2 - Ajtai Commitment

## Step 3 - Proyections (Johnson-Lindenstrauss) 

## Step 4 - Aggregation

## Step 5 - Amortization  

