pub fn say_hello() {
    println!("Hello, world!");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_say_hello() {
        say_hello();
    }

    #[test]
    fn test_for_workflows() {
        my_func(&"world"); // Should be a Clippy warning
        assert_eq!(true, false)
    }

    fn my_func(s: &str) {
        println!("Hello, {s}!");
    }
}
