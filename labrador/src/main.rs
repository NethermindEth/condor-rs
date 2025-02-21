use labrador::zq::Zq;

fn main() {
    let a = Zq::new(5);
    let b = Zq::new(3);
    println!("a + b = {}", a + b);
}
