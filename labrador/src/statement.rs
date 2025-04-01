use crate::poly::PolyVector;

/// Statement is the input of the prover, which contains the constraints and the witness.
/// All parameters are from line 1, st(), in the verifier process, page 18 from the paper.
pub struct Statement {
    // $a_{ij}^k$
    pub a_constraint: Vec<Vec<PolyVector>>,
    // $\varphi_i^k$
    pub phi_constraint: Vec<Vec<PolyVector>>,
    // $b^{(k)}$
    pub b_constraint: PolyVector,
    // $a_{ij}^{'(l)}$
    pub a_ct: Vec<Vec<PolyVector>>,
    // $\varphi_i^{'(l)}$
    pub phi_ct: Vec<Vec<PolyVector>>,
    // $b_0^{'(l)}$
    pub b_ct: PolyVector,
}
