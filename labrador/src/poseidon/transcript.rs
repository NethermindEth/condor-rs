use crate::zq::Zq;

use super::{permutation::PoseidonPermutation, sponge::PoseidonSponge};

// ===========================================================================
// Poseidon parameter set – hard-coded for Fiat–Shamir transcripts
// ===========================================================================
//
// The ARK (add-round-key) and MDS matrices below were *generated once* and are
// embedded so that every prover/verifier derives the **same** challenges.
const WIDTH: usize = 9;
const RATE: usize = 1;
const OUTPUT_LEN: usize = 1;
const ROUNDS: usize = 28;
const PARTIAL_ROUNDS: usize = 20;
const ALPHA: u64 = 3;
// Hard‑coded round constants (ARK) and MDS matrix.  Generated off‑line.
const ARK: [[Zq; WIDTH]; ROUNDS] = [
    [
        Zq::new(2773654128),
        Zq::new(2698871901),
        Zq::new(1186535217),
        Zq::new(2276367172),
        Zq::new(460986456),
        Zq::new(1130506256),
        Zq::new(1436990685),
        Zq::new(3571618093),
        Zq::new(1928846831),
    ],
    [
        Zq::new(784272210),
        Zq::new(2362493121),
        Zq::new(3754859051),
        Zq::new(1152838762),
        Zq::new(1277027414),
        Zq::new(4290089713),
        Zq::new(750711681),
        Zq::new(3735304872),
        Zq::new(2484390012),
    ],
    [
        Zq::new(3022291695),
        Zq::new(1296513440),
        Zq::new(2666656071),
        Zq::new(1392772787),
        Zq::new(3924839390),
        Zq::new(3434015439),
        Zq::new(3547712534),
        Zq::new(3119852137),
        Zq::new(3798005745),
    ],
    [
        Zq::new(2375363973),
        Zq::new(1758886039),
        Zq::new(2041597959),
        Zq::new(3263939143),
        Zq::new(1885616341),
        Zq::new(1813930986),
        Zq::new(203888102),
        Zq::new(1528124030),
        Zq::new(56779178),
    ],
    [
        Zq::new(3755261782),
        Zq::new(706143006),
        Zq::new(2082983570),
        Zq::new(352354306),
        Zq::new(2353150435),
        Zq::new(3650061737),
        Zq::new(3893118498),
        Zq::new(183150537),
        Zq::new(2228089161),
    ],
    [
        Zq::new(922618003),
        Zq::new(3292151780),
        Zq::new(236167017),
        Zq::new(2617694273),
        Zq::new(2876369390),
        Zq::new(4265817939),
        Zq::new(1383107438),
        Zq::new(286389424),
        Zq::new(3869373395),
    ],
    [
        Zq::new(1955724147),
        Zq::new(1111197896),
        Zq::new(633124926),
        Zq::new(2523587424),
        Zq::new(4135563482),
        Zq::new(3059457948),
        Zq::new(2282176497),
        Zq::new(1783995730),
        Zq::new(3403269348),
    ],
    [
        Zq::new(508988721),
        Zq::new(301251615),
        Zq::new(4208429837),
        Zq::new(2847939307),
        Zq::new(3727165119),
        Zq::new(3625918486),
        Zq::new(3478562333),
        Zq::new(1459737698),
        Zq::new(100323470),
    ],
    [
        Zq::new(3179318671),
        Zq::new(1533355377),
        Zq::new(2704441550),
        Zq::new(2362155401),
        Zq::new(4225681297),
        Zq::new(621873557),
        Zq::new(1849485098),
        Zq::new(4012360807),
        Zq::new(3322683455),
    ],
    [
        Zq::new(271837745),
        Zq::new(1175116485),
        Zq::new(1502825906),
        Zq::new(113799952),
        Zq::new(710728204),
        Zq::new(930285870),
        Zq::new(936814998),
        Zq::new(3142336897),
        Zq::new(1195857553),
    ],
    [
        Zq::new(3229477597),
        Zq::new(3874900104),
        Zq::new(633468355),
        Zq::new(664022320),
        Zq::new(1037768339),
        Zq::new(24745787),
        Zq::new(1390286585),
        Zq::new(2050114210),
        Zq::new(3461313472),
    ],
    [
        Zq::new(2038116502),
        Zq::new(2756126285),
        Zq::new(314702156),
        Zq::new(198049006),
        Zq::new(3362107641),
        Zq::new(804124193),
        Zq::new(1347345692),
        Zq::new(3778418349),
        Zq::new(1603517811),
    ],
    [
        Zq::new(246075735),
        Zq::new(2910920140),
        Zq::new(2101057122),
        Zq::new(2041317133),
        Zq::new(2883309822),
        Zq::new(4019460306),
        Zq::new(2468923228),
        Zq::new(964910736),
        Zq::new(1333043724),
    ],
    [
        Zq::new(293190938),
        Zq::new(3902017027),
        Zq::new(549889505),
        Zq::new(4272274465),
        Zq::new(3832846664),
        Zq::new(2204806555),
        Zq::new(3020552376),
        Zq::new(3447507883),
        Zq::new(855953368),
    ],
    [
        Zq::new(2079971112),
        Zq::new(2508401160),
        Zq::new(3834351883),
        Zq::new(1427569236),
        Zq::new(3358408800),
        Zq::new(2439693552),
        Zq::new(4054927129),
        Zq::new(4053693217),
        Zq::new(723473858),
    ],
    [
        Zq::new(1025339977),
        Zq::new(585257694),
        Zq::new(3070506491),
        Zq::new(3175819591),
        Zq::new(1392016786),
        Zq::new(1516602977),
        Zq::new(3228410444),
        Zq::new(3130776340),
        Zq::new(3165959881),
    ],
    [
        Zq::new(3842311420),
        Zq::new(2262157250),
        Zq::new(3756225800),
        Zq::new(2582280172),
        Zq::new(1270243441),
        Zq::new(2052508852),
        Zq::new(4271664205),
        Zq::new(620984828),
        Zq::new(3325253760),
    ],
    [
        Zq::new(548636457),
        Zq::new(1164607978),
        Zq::new(3027080685),
        Zq::new(1908414320),
        Zq::new(1208816014),
        Zq::new(671096184),
        Zq::new(4042441851),
        Zq::new(2000153875),
        Zq::new(2975622655),
    ],
    [
        Zq::new(2588519128),
        Zq::new(1500421688),
        Zq::new(280920169),
        Zq::new(1425318099),
        Zq::new(173722839),
        Zq::new(1354484516),
        Zq::new(1969077869),
        Zq::new(2165032598),
        Zq::new(2764691127),
    ],
    [
        Zq::new(3889006152),
        Zq::new(676568773),
        Zq::new(1097350839),
        Zq::new(3748044660),
        Zq::new(3451456889),
        Zq::new(1469067299),
        Zq::new(1196266786),
        Zq::new(3351660008),
        Zq::new(3216956496),
    ],
    [
        Zq::new(2920308927),
        Zq::new(1964475343),
        Zq::new(1208178427),
        Zq::new(1685596759),
        Zq::new(3452869029),
        Zq::new(574413953),
        Zq::new(1372297329),
        Zq::new(1642685791),
        Zq::new(1521032998),
    ],
    [
        Zq::new(1396137219),
        Zq::new(3249324467),
        Zq::new(1491626859),
        Zq::new(1222034466),
        Zq::new(3741121943),
        Zq::new(2630086446),
        Zq::new(4101711170),
        Zq::new(3962365070),
        Zq::new(630331971),
    ],
    [
        Zq::new(3756693713),
        Zq::new(1976849786),
        Zq::new(1599868415),
        Zq::new(4193484663),
        Zq::new(3575428569),
        Zq::new(2614894763),
        Zq::new(721681014),
        Zq::new(3645252436),
        Zq::new(96433812),
    ],
    [
        Zq::new(1889643081),
        Zq::new(178898262),
        Zq::new(2836025067),
        Zq::new(1086358336),
        Zq::new(3421342207),
        Zq::new(1160658413),
        Zq::new(3078690548),
        Zq::new(1734238039),
        Zq::new(1684918153),
    ],
    [
        Zq::new(4186013047),
        Zq::new(1024422138),
        Zq::new(4025507495),
        Zq::new(2413389692),
        Zq::new(614405915),
        Zq::new(2631560766),
        Zq::new(4144324857),
        Zq::new(2400759460),
        Zq::new(1279501883),
    ],
    [
        Zq::new(1791688840),
        Zq::new(2006611136),
        Zq::new(3093261711),
        Zq::new(1016157743),
        Zq::new(11561707),
        Zq::new(1430464813),
        Zq::new(73038415),
        Zq::new(2025372355),
        Zq::new(2402576570),
    ],
    [
        Zq::new(2758798058),
        Zq::new(2132722282),
        Zq::new(1698416471),
        Zq::new(3578176273),
        Zq::new(550672905),
        Zq::new(70724481),
        Zq::new(3070901761),
        Zq::new(306186726),
        Zq::new(2828652630),
    ],
    [
        Zq::new(1005410618),
        Zq::new(1376511878),
        Zq::new(4093260411),
        Zq::new(3950256285),
        Zq::new(1890541837),
        Zq::new(1315511419),
        Zq::new(3327986955),
        Zq::new(1322306885),
        Zq::new(1270447424),
    ],
];
const MDS: [[Zq; WIDTH]; WIDTH] = [
    [
        Zq::new(3997651091),
        Zq::new(658253063),
        Zq::new(2095400994),
        Zq::new(440926105),
        Zq::new(1796368741),
        Zq::new(1961349520),
        Zq::new(540996892),
        Zq::new(35935778),
        Zq::new(732075401),
    ],
    [
        Zq::new(1743127377),
        Zq::new(3247221626),
        Zq::new(4148659237),
        Zq::new(2205285832),
        Zq::new(3910687740),
        Zq::new(2853895742),
        Zq::new(1662369266),
        Zq::new(241107308),
        Zq::new(4022943497),
    ],
    [
        Zq::new(2071371317),
        Zq::new(4222913952),
        Zq::new(1163051760),
        Zq::new(3185249567),
        Zq::new(4114377118),
        Zq::new(1048229747),
        Zq::new(658207529),
        Zq::new(1971020081),
        Zq::new(4001507915),
    ],
    [
        Zq::new(4278306557),
        Zq::new(2902580578),
        Zq::new(1382067936),
        Zq::new(3738299059),
        Zq::new(1548047067),
        Zq::new(1128497469),
        Zq::new(2358825001),
        Zq::new(567535302),
        Zq::new(4263023986),
    ],
    [
        Zq::new(188724725),
        Zq::new(890882199),
        Zq::new(2995073197),
        Zq::new(3520401121),
        Zq::new(2334136676),
        Zq::new(1039560330),
        Zq::new(3354448203),
        Zq::new(2456760197),
        Zq::new(860596610),
    ],
    [
        Zq::new(1063803022),
        Zq::new(4115470744),
        Zq::new(2914526459),
        Zq::new(995387934),
        Zq::new(47066002),
        Zq::new(3823938923),
        Zq::new(329509446),
        Zq::new(2091055303),
        Zq::new(2927144621),
    ],
    [
        Zq::new(1222968641),
        Zq::new(3572964758),
        Zq::new(532234855),
        Zq::new(1137915407),
        Zq::new(902452800),
        Zq::new(1725021750),
        Zq::new(3941355778),
        Zq::new(816646153),
        Zq::new(2658636777),
    ],
    [
        Zq::new(2652215134),
        Zq::new(3067867604),
        Zq::new(1764433633),
        Zq::new(3289297051),
        Zq::new(3187448683),
        Zq::new(3735943544),
        Zq::new(3993652112),
        Zq::new(2905114558),
        Zq::new(477157260),
    ],
    [
        Zq::new(3093659876),
        Zq::new(479575061),
        Zq::new(3390824282),
        Zq::new(3037823443),
        Zq::new(1492293764),
        Zq::new(3797464643),
        Zq::new(3579348020),
        Zq::new(2374361921),
        Zq::new(842530663),
    ],
];

/// Convenience type aliases for the fixed permutation + sponge
type Perm = PoseidonPermutation<WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>;
type Sponge = PoseidonSponge<OUTPUT_LEN, RATE, WIDTH, ROUNDS, PARTIAL_ROUNDS, ALPHA>;

/// A **Fiat–Shamir transcript** based on Poseidon.
///
/// The transcript collects field elements via [`absorb_element`].  When a new
/// challenge is required, [`get_scalar_challenge`] returns the Poseidon hash of
/// the current buffer, interpreted as a scalar in the same field `Zq`.
///
/// The internal message buffer is **not** cleared after squeezing, allowing the
/// caller to derive *multiple* independent challenges sequentially by first
/// appending the previous challenge back into the transcript (standard
/// practice in many proof systems).
pub struct Transcript {
    buf: Vec<Zq>,
}

impl Default for Transcript {
    fn default() -> Self {
        Self::new()
    }
}

impl Transcript {
    /// Creates an **empty** transcript.
    pub fn new() -> Self {
        Self { buf: Vec::new() }
    }

    /// Appends a field element to the transcript.  Elements are stored in the
    /// order they are seen and fed verbatim to the Poseidon sponge.
    pub fn absorb_element(&mut self, elem: Zq) {
        self.buf.push(elem);
    }

    /// Generates a *scalar* challenge and **appends it** to the transcript so
    /// that subsequent challenges are bound to the previous ones.
    pub fn get_scalar_challenge(&mut self) -> Zq {
        let permutation = Perm::new_with_ark_mds(ARK, MDS);
        let mut sponge = Sponge::new(self.buf.clone(), permutation);
        let challenge = sponge.compute_hash()[0];
        self.buf.push(challenge);
        challenge
    }
}

#[cfg(test)]
mod tests {
    use super::Transcript;
    use crate::zq::Zq;

    #[test]
    fn test_transcript() {
        let mut transcript = Transcript::new();
        transcript.absorb_element(Zq::new(12));
        transcript.absorb_element(Zq::new(12312));
        transcript.absorb_element(Zq::new(111));
        let _result = transcript.get_scalar_challenge();
        // result == 3770429813
    }

    #[test]
    fn test_challenge_is_appended_to_transcript() {
        let mut t = Transcript::new();
        t.absorb_element(Zq::new(42));
        let c1 = t.get_scalar_challenge();
        let c2 = t.get_scalar_challenge();
        assert_ne!(c1, c2);
    }

    /// Two transcripts with identical message sequences must agree on the
    /// challenge value.
    #[test]
    fn test_two_transcripts_agree() {
        let mut t1 = Transcript::new();
        let mut t2 = Transcript::new();
        for &x in &[1u32, 2, 3, 5, 8] {
            t1.absorb_element(Zq::new(x));
            t2.absorb_element(Zq::new(x));
        }
        assert_eq!(t1.get_scalar_challenge(), t2.get_scalar_challenge());
    }

    /// Changing *one* absorbed element should almost certainly alter the
    /// challenge (collision resistance is probabilistic; here we just assert
    /// inequality for this specific test vector).
    #[test]
    fn test_challenge_changes_with_input() {
        let mut t1 = Transcript::new();
        let mut t2 = Transcript::new();
        t1.absorb_element(Zq::new(99));
        t2.absorb_element(Zq::new(100));
        assert_ne!(t1.get_scalar_challenge(), t2.get_scalar_challenge());
    }
}
