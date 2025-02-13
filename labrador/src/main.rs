use labrador::say_hello;
use labrador::zq::Zq;

mod poly_ring;
use poly_ring::Poly;

fn main() {
    say_hello();

    // Example poly_ring
    const D: usize = 3; // Define the constant d in (x^d + 1)
    let p1 = Poly::<D>::create_poly(vec![1, 2]);
    let p2 = Poly::<D>::create_poly(vec![1]);

    // Multiply polynomials
    let product = p1.mul(&p2);
    // Add polynomials
    let sum = p1.add(&p2);
    // Dot product between coefficients
    let dot = p1.inner_product(&p2);
    println!("Product: {:?}", product);
    println!("sum: {:?}", sum);
    println!("dot: {:?}", dot);

    let a = Zq::new(5);
    let b = Zq::new(3);
    println!("a + b = {}", a + b);

}
