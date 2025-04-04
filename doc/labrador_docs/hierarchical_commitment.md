# Hierarchical Commitment Structure

A key optimization in the LaBRADOR protocol is the use of hierarchical commitments - a two-layer structure with inner and outer commitments.

## Inner and Outer Commitments

The hierarchical commitment structure consists of:

1. **Inner Commitments**: First-level commitments to witness vectors
2. **Outer Commitments**: Second-level commitments that summarize multiple inner commitments

This two-layer structure is crucial for reducing proof size by allowing the prover to move material from the current proof iteration to the next iteration, where it will be further compressed.

## Decomposition Process

When creating a hierarchical commitment:

1. The prover first commits to each witness vector using individual Ajtai commitments, creating the inner commitments
2. These inner commitments are then decomposed into a base-B representation to reduce coefficient size
3. The decomposed parts are concatenated into a single vector
4. This vector is then committed to using another Ajtai commitment, creating the outer commitment

The decomposition uses centered representatives modulo the base to minimize the magnitude of coefficients, ensuring compactness.

## Verification Process

To verify a hierarchical commitment:

1. The verifier checks each inner commitment against its opening information
2. The verifier then decomposes the inner commitments using the same process
3. The decomposed parts are concatenated and committed to using the outer commitment scheme
4. The resulting commitment is compared against the provided outer commitment
5. The outer commitment's opening information is also verified

## Benefits of Hierarchical Commitments

The hierarchical structure provides several advantages:

1. **Reduced Communication Complexity**: By summarizing multiple inner commitments into a single outer commitment, the total size of the proof is reduced
2. **Efficient Verification**: The verification process can be performed efficiently even with the two-layer structure
3. **Flexibility**: The decomposition parameters can be adjusted to optimize for specific use cases

## Integration with LaBRADOR Protocol

In the context of the LaBRADOR protocol, hierarchical commitments allow moving material from π(i+1) to s(i+1), which is then shrunk in subsequent iterations, resulting in a more compact overall proof size.