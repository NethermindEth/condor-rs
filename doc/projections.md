## Step 3 - Projections (Johnson-Lindenstrauss)
The main idea behind the Modular Johnson-Lindenstrauss Lemma is that if you want to prove that a long vector has a small norm, instead of having to send the entire vector to verify this condition, the verifier can send random linear projections as challenges to the prover and then the prover just needs to return the aplied proyection to the vector. The lemma guarantees that the projections almost preserve the L2 norm, so by performing a sufficient number of challenges, one can gain information about the actual norm size of the original vector.

