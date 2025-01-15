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
    fn test_for_workflow() {
        assert_eq!(true, false)
    }
}
