pub trait Transcript {
    fn new() -> Self;
    fn absorb(&mut self, value: u64);
    fn get_challenge(&self) -> u64;
}
pub struct SimpleTranscript {
    state: u64,
}
impl Transcript for SimpleTranscript {
    fn new() -> Self {
        Self { state: 0 }  // Start with 0
    }

    fn absorb(&mut self, value: u64) {
        self.state ^= value;  // Just XOR for simplicity
    }

    fn get_challenge(&self) -> u64 {
        self.state.wrapping_add(1)  // Return state + 1 as a dummy challenge
    }
}