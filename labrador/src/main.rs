// use labrador::core::jl::{verify_upper_bound, ProjectionMatrix, ProjectionVector};
// use labrador::ring::rq::Rq;
// use labrador::ring::rq_vector::RqVector;
// use labrador::ring::zq::Zq;
// use rand::rng;

// const D: usize = 4; // Degree of polynomials in S_i
// const N: usize = 5; // Size of S_i

// fn main() {
//     // Example poly_ring
//     let p1: Rq<D> = vec![Zq::new(1)].into();
//     let p2: Rq<D> = vec![Zq::new(2), Zq::new(1), Zq::new(1)].into();
//     // Perform polynomial multiplication
//     let product = p1.clone() * p2.clone();

//     // Perform polynomial addition
//     let sum = p1.clone() + p2.clone();

//     // Perform polynomial subtraction
//     let sub = p1.clone() - p2.clone();

//     // Compute the dot product between the polynomial coefficients
//     let dot = p1.clone().inner_product(&p2);

//     // Negate the polynomial
//     let negation = -p1.clone();

//     // Perform scalar multiplication
//     let scalar_multiplication = p1.scalar_mul(Zq::new(2));

//     // Evaluate the polynomial at x = 2
//     let evaluation = p2.eval(Zq::new(2));

//     // Check if the polynomial is the zero polynomial
//     let zero_check = p1.is_zero();

//     // Check if p1 is equal to p2
//     let are_equal = p1.is_equal(&p2);

//     // Print the results
//     println!("Product: {:?}", product);
//     println!("Sum: {:?}", sum);
//     println!("Subtraction: {:?}", sub);
//     println!("Dot product: {:?}", dot);
//     println!("Negation: {:?}", negation);
//     println!("Scalar multiplication: {:?}", scalar_multiplication);
//     println!("Evaluation at x=2: {:?}", evaluation);
//     println!("Is zero polynomial: {:?}", zero_check);
//     println!("Are polynomials equal: {:?}", are_equal);

//     let a = Zq::new(5);
//     let b = Zq::new(3);
//     println!("a + b = {}", a + b);

//     // Johnson Linderstrauss Projections
//     // Example
//     // Generate the random polynomials
//     let n = 3;
//     let mut rng = rng();
//     let polynomials = RqVector::<N, D>::random_ternary(&mut rng);
//     // Random projection matrix
//     let matrix = ProjectionMatrix::new(n);
//     // Calculate projection
//     let projection = ProjectionVector::new(&matrix, &polynomials);
//     // Within bounds with probability 1/2
//     let beta = polynomials.compute_norm_squared();
//     println!("{}", verify_upper_bound(projection, beta));
// }

fn main() {}
