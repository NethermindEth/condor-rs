use labrador::rq::Rq;
use labrador::say_hello;
use labrador::zq::Zq;

fn main() {
    say_hello();

    // Example poly_ring
    let p1: Rq<2> = vec![Zq::new(1)].into();
    let p2: Rq<2> = vec![Zq::new(2), Zq::new(1), Zq::new(1)].into();
    // Perform polynomial multiplication
    let product = p1 * p2;

    // Perform polynomial addition
    let sum = p1 + p2;

    // Perform polynomial subtraction
    let sub = p1 - p2;

    // Compute the dot product between the polynomial coefficients
    let dot = p1.inner_product(&p2);

    // Negate the polynomial
    let negation = p1.neg();

    // Perform scalar multiplication
    let scalar_multiplication = p1.scalar_mul(Zq::new(2));

    // Evaluate the polynomial at x = 2
    let evaluation = p1.eval(Zq::new(2));

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
}
