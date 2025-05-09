use crate::ring::rq::Rq;
use super::shake_sponge::Sponge;

pub struct LabradorTranscript<S: Sponge> {
    sponge: S,
}

impl<S: Sponge> LabradorTranscript<S> {
    fn new(sponge: S) -> Self {
        Self {
            sponge,
        }
    }

    fn absorb(&mut self, ring: &[Rq]) {
        self.sponge.absorb(ring);
    }

    fn get_challenges(&mut self, length: usize) -> Vec<Rq> {
        let challenge = self.sponge.squeeze(length);
        challenge
    }
}

#[cfg(test)]
mod tests {
    use crate::{poseidon::shake_sponge::ShakeSponge, ring::zq::Zq};

    use super::*;

    #[test]
    fn test_transcript_is_deterministic() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);
        let polynomial_3 = Rq::new([Zq::new(9891741); Rq::DEGREE]);

        let mut transcript1 = LabradorTranscript::new(ShakeSponge::default());
        transcript1.absorb(&[polynomial_1, polynomial_2, polynomial_3]);

        let mut transcript2 = LabradorTranscript::new(ShakeSponge::default());
        transcript2.absorb(&[polynomial_1, polynomial_2, polynomial_3]);
        assert_eq!(transcript1.get_challenges(1), transcript2.get_challenges(1))
    }

    #[test]
    fn test_transcript_execution() {
        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        let polynomial_1 = Rq::new([Zq::new(2); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(5); Rq::DEGREE]);
        let polynomial_3 = Rq::new([Zq::new(1); Rq::DEGREE]);

        transcript.absorb(&[polynomial_1, polynomial_2, polynomial_3]);
        let result = transcript.get_challenges(1);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_transcript_message_order() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);

        let mut transcript1 = LabradorTranscript::new(ShakeSponge::default());
        transcript1.absorb(&[polynomial_1, polynomial_2]);

        let mut transcript2 = LabradorTranscript::new(ShakeSponge::default());
        transcript2.absorb(&[polynomial_2, polynomial_1]);
        assert_ne!(transcript1.get_challenges(1), transcript2.get_challenges(1))
    }

    #[test]
    fn test_challenge_is_appended() {
        let polynomial_1 = Rq::new([Zq::new(54821); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(2131213); Rq::DEGREE]);

        let mut transcript = LabradorTranscript::new(ShakeSponge::default());
        transcript.absorb(&[polynomial_1, polynomial_2]);
        let result1 = transcript.get_challenges(1);
        let result2 = transcript.get_challenges(1);
        assert_ne!(result1, result2)
    }

    #[test]
    fn test_output_length_correctness() {
        let sponge = ShakeSponge::default();
        let mut transcript = LabradorTranscript::new(sponge);
        let polynomial_1 = Rq::new([Zq::new(2); Rq::DEGREE]);
        let polynomial_2 = Rq::new([Zq::new(5); Rq::DEGREE]);
        let polynomial_3 = Rq::new([Zq::new(1); Rq::DEGREE]);
        let polynomial_4 = Rq::new([Zq::new(9); Rq::DEGREE]);
        transcript.absorb(&[polynomial_1, polynomial_2, polynomial_3, polynomial_4]);

        let result = transcript.get_challenges(2);
        assert_eq!(result.len(), 2);
    }
}