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
        use base64::prelude::*;
        let _ = BASE64_STANDARD.encode([1,2,3,4,5]);
        assert_eq!(true, true)
    }
}
