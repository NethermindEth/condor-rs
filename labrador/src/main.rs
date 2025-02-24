use labrador::jl::compute_norm;
use labrador::jl::generate_random_polynomials;
use labrador::jl::ProjectionMatrix;
use labrador::jl::ProjectionVector;
use labrador::rq::Rq;
use labrador::zq::Zq;
use rand::rngs::ThreadRng; // Import ThreadRng

const D: usize = 4; // Degree of polynomials (you can change this to the required degree)
const N: usize = 3; // Number of polynomials (you can adjust this as needed)

fn main() {
    // Example poly_ring
    let p1: Rq<2> = vec![Zq::new(1)].into();
    let p2: Rq<2> = vec![Zq::new(2), Zq::new(1), Zq::new(1)].into();
    // Perform polynomial multiplication
    let product = p1.clone() * p2.clone();

    // Perform polynomial addition
    let sum = p1.clone() + p2.clone();

    // Perform polynomial subtraction
    let sub = p1.clone() - p2.clone();

    // Compute the dot product between the polynomial coefficients
    let dot = p1.clone().inner_product(&p2);

    // Negate the polynomial
    let negation = -p1.clone();

    // Perform scalar multiplication
    let scalar_multiplication = p1.scalar_mul(Zq::new(2));

    // Evaluate the polynomial at x = 2
    let evaluation = p2.eval(Zq::new(2));

    // Check if the polynomial is the zero polynomial
    let zero_check = p1.is_zero();

    // Check if p1 is equal to p2
    let are_equal = p1.is_equal(&p2);

    // Print the results
    println!("Product: {:?}", product);
    println!("Sum: {:?}", sum);
    println!("Subtraction: {:?}", sub);
    println!("Dot product: {:?}", dot);
    println!("Negation: {:?}", negation);
    println!("Scalar multiplication: {:?}", scalar_multiplication);
    println!("Evaluation at x=2: {:?}", evaluation);
    println!("Is zero polynomial: {:?}", zero_check);
    println!("Are polynomials equal: {:?}", are_equal);

    let a = Zq::new(5);
    let b = Zq::new(3);
    println!("a + b = {}", a + b);

    // Johnson Linderstrauss Projections
    // Example parameters

    // Generate the random polynomials
    let beta: Zq = Zq::new(500);
    let mut rng = rand::rngs::ThreadRng::default();
    let polynomials = generate_random_polynomials::<ThreadRng, D>(N, &mut rng, beta);
    let matrix = ProjectionMatrix::new(D * N);
    let projection = ProjectionVector::new(&matrix, &polynomials);

    // Print the generated polynomial Norms
    println!(
        "beta = {} | Polynomial Norm = {} | sqrt(128) * norm polynomials = {} | Projection Norm  = {}",beta,
        compute_norm(&polynomials),Zq::new(128) * compute_norm(&polynomials),projection.norm()
    );

    // Show Polynomials S_i
    //for (i, poly) in polynomials.iter().enumerate() {
    //    println!("polynomial {} = {:?}", i, poly);
    //}

    // Show projection elements
    //println!("elements = {:?}", projection.get_projection())
}
