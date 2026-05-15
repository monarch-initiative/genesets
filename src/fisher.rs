#[derive(Clone, Debug)]
pub struct FisherCache {
    log_factorials: Vec<f64>,
}

impl FisherCache {
    pub fn new(max_n: usize) -> Self {
        let mut log_factorials = Vec::with_capacity(max_n + 1);
        log_factorials.push(0.0);
        for n in 1..=max_n {
            let next = log_factorials[n - 1] + (n as f64).ln();
            log_factorials.push(next);
        }
        Self { log_factorials }
    }

    pub fn right_tail(&self, population: u32, successes: u32, draws: u32, observed: u32) -> f64 {
        if population == 0 || draws > population || successes > population {
            return f64::NAN;
        }
        if observed == 0 {
            return 1.0;
        }

        let max_observed = successes.min(draws);
        let min_observed = draws.saturating_sub(population - successes);
        if observed <= min_observed {
            return 1.0;
        }
        if observed > max_observed {
            return 0.0;
        }

        let log_p = self.log_hypergeom(population, successes, draws, observed);
        let mut p = log_p.exp();
        if p == 0.0 {
            return self.right_tail_log_sum(population, successes, draws, observed, max_observed);
        }

        let mut sum = p;
        for i in observed..max_observed {
            let numerator = (successes - i) as f64 * (draws - i) as f64;
            let second_denominator =
                population as i64 - successes as i64 - draws as i64 + i as i64 + 1;
            let denominator = (i + 1) as f64 * second_denominator as f64;
            if numerator == 0.0 || denominator <= 0.0 {
                break;
            }
            p *= numerator / denominator;
            sum += p;
        }
        sum.clamp(0.0, 1.0)
    }

    fn right_tail_log_sum(
        &self,
        population: u32,
        successes: u32,
        draws: u32,
        observed: u32,
        max_observed: u32,
    ) -> f64 {
        let logs: Vec<f64> = (observed..=max_observed)
            .map(|i| self.log_hypergeom(population, successes, draws, i))
            .collect();
        let max_log = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !max_log.is_finite() {
            return 0.0;
        }
        let scaled_sum: f64 = logs.iter().map(|log| (log - max_log).exp()).sum();
        (max_log + scaled_sum.ln()).exp().clamp(0.0, 1.0)
    }

    fn log_hypergeom(&self, population: u32, successes: u32, draws: u32, observed: u32) -> f64 {
        self.log_choose(successes, observed)
            + self.log_choose(population - successes, draws - observed)
            - self.log_choose(population, draws)
    }

    fn log_choose(&self, n: u32, k: u32) -> f64 {
        if k > n {
            return f64::NEG_INFINITY;
        }
        let n = n as usize;
        let k = k.min(n as u32 - k) as usize;
        self.log_factorials[n] - self.log_factorials[k] - self.log_factorials[n - k]
    }
}

#[cfg(test)]
mod tests {
    use super::FisherCache;

    #[test]
    fn computes_known_right_tail() {
        let fisher = FisherCache::new(5);
        let p = fisher.right_tail(5, 2, 3, 2);
        assert!((p - 0.3).abs() < 1e-12);

        let p = fisher.right_tail(5, 3, 3, 2);
        assert!((p - 0.7).abs() < 1e-12);
    }

    #[test]
    fn handles_impossible_or_trivial_tails() {
        let fisher = FisherCache::new(10);
        assert_eq!(fisher.right_tail(10, 3, 4, 0), 1.0);
        assert_eq!(fisher.right_tail(10, 3, 4, 5), 0.0);
    }
}
