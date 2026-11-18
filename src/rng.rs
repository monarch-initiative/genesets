//! Small reproducible random number generator shared by the samplers.

/// xoshiro256++ seeded through SplitMix64; small, fast, and reproducible.
pub(crate) struct Rng {
    state: [u64; 4],
}

impl Rng {
    pub(crate) fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Self {
            state: [next(), next(), next(), next()],
        }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        let s = &mut self.state;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    pub(crate) fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    pub(crate) fn index(&mut self, len: usize) -> usize {
        ((u128::from(self.next_u64()) * len as u128) >> 64) as usize
    }

    pub(crate) fn weighted_index(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut remaining = self.uniform() * total;
        for (index, weight) in weights.iter().enumerate() {
            if *weight > 0.0 && remaining < *weight {
                return index;
            }
            remaining -= weight;
        }
        weights
            .iter()
            .rposition(|weight| *weight > 0.0)
            .unwrap_or(0)
    }
}
