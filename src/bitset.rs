#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DenseBitSet {
    words: Vec<u64>,
    len: usize,
}

impl DenseBitSet {
    pub fn new(len: usize) -> Self {
        let word_count = len.div_ceil(64);
        Self {
            words: vec![0; word_count],
            len,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn set(&mut self, index: usize) {
        debug_assert!(index < self.len);
        let word = index / 64;
        let bit = index % 64;
        self.words[word] |= 1u64 << bit;
    }

    pub fn contains(&self, index: usize) -> bool {
        debug_assert!(index < self.len);
        let word = index / 64;
        let bit = index % 64;
        (self.words[word] & (1u64 << bit)) != 0
    }

    pub fn count_ones(&self) -> u32 {
        self.words.iter().map(|word| word.count_ones()).sum()
    }

    pub fn or_with(&mut self, other: &Self) {
        self.assert_same_len(other);
        for (left, right) in self.words.iter_mut().zip(&other.words) {
            *left |= right;
        }
    }

    pub fn and_count(&self, other: &Self) -> u32 {
        self.assert_same_len(other);
        self.words
            .iter()
            .zip(&other.words)
            .map(|(left, right)| (*left & *right).count_ones())
            .sum()
    }

    pub fn and3_count(&self, other: &Self, third: &Self) -> u32 {
        self.assert_same_len(other);
        self.assert_same_len(third);
        self.words
            .iter()
            .zip(&other.words)
            .zip(&third.words)
            .map(|((left, right), third)| (*left & *right & *third).count_ones())
            .sum()
    }

    pub fn intersection_indices3(&self, other: &Self, third: &Self) -> Vec<usize> {
        self.assert_same_len(other);
        self.assert_same_len(third);

        let mut indices = Vec::new();
        for (word_index, ((left, right), third)) in self
            .words
            .iter()
            .zip(&other.words)
            .zip(&third.words)
            .enumerate()
        {
            let mut word = *left & *right & *third;
            while word != 0 {
                let bit = word.trailing_zeros() as usize;
                let index = word_index * 64 + bit;
                if index < self.len {
                    indices.push(index);
                }
                word &= word - 1;
            }
        }
        indices
    }

    fn assert_same_len(&self, other: &Self) {
        debug_assert_eq!(self.len, other.len);
        debug_assert_eq!(self.words.len(), other.words.len());
    }
}

#[cfg(test)]
mod tests {
    use super::DenseBitSet;

    #[test]
    fn counts_intersections() {
        let mut left = DenseBitSet::new(130);
        let mut right = DenseBitSet::new(130);
        let mut bg = DenseBitSet::new(130);

        for index in [1, 64, 65, 129] {
            left.set(index);
        }
        for index in [0, 64, 65, 80, 129] {
            right.set(index);
        }
        for index in [64, 65, 80] {
            bg.set(index);
        }

        assert_eq!(left.count_ones(), 4);
        assert_eq!(left.and_count(&right), 3);
        assert_eq!(left.and3_count(&right, &bg), 2);
        assert_eq!(left.intersection_indices3(&right, &bg), vec![64, 65]);
    }
}
