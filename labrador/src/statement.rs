use crate::poly::PolyVector;

/// Statement is the input of the prover, which contains the constraints and the witness.
/// All parameters are from equation 1, in the verifier process, page 18 from the paper.
pub struct Statement {
    pub a_constraint: Vec<Vec<PolyVector>>,
    pub phi_constraint: Vec<Vec<PolyVector>>,
    pub b_constraint: PolyVector,
    pub a_ct: Vec<Vec<PolyVector>>,
    pub phi_ct: Vec<Vec<PolyVector>>,
    pub b_ct: PolyVector,
}
