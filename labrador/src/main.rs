use labrador::say_hello;
use labrador::zq::Zq;

fn main() {
    say_hello();
    let a = Zq::new(5);
    let b = Zq::new(3);
    println!("a + b = {}", a + b);
}
