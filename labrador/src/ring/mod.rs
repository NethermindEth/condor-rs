pub mod rq;
pub mod rq_matrix;
pub mod rq_vector;
pub mod zq;

pub trait Norms {
    type NormType;
    fn l2_norm_squared(&self) -> Self::NormType;
    fn linf_norm(&self) -> Self::NormType;
}
