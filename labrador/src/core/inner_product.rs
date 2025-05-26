use std::{
    borrow::Borrow,
    ops::{Add, Mul},
};

pub fn compute_linear_combination<E, B, C>(elements: &[E], challenges: &[C]) -> B
where
    E: Borrow<B>,
    for<'a> &'a B: Mul<&'a C, Output = B>,
    for<'a> &'a B: Add<&'a B, Output = B>,
{
    debug_assert_eq!(
        elements.len(),
        challenges.len(),
        "vectors must be the same length"
    );
    debug_assert!(!elements.is_empty(), "`elements` must not be empty");

    let mut zipped_iter = elements.iter().zip(challenges.iter());
    // Must do the following as the init value in fold requires size of B
    let (e0, c0) = zipped_iter.next().unwrap();
    let init = e0.borrow() * c0;

    zipped_iter.fold(init, |acc, (elem, c)| &acc + &(elem.borrow() * c))
}
