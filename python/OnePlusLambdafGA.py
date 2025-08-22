

import numpy as np
import ioh

class OnePlusLambdafGA:
    """
    (1+λ) Evolutionary Algorithm for bitstring domains (e.g., IOHprofiler PBO).

    Usage:
        algo = OnePlusLambdaEA(lam=10)
        fx_best, x_best = algo(problem)   # where `problem` is an IOH function object

    Behavior:
      - Maintains a single parent x
      - Each iteration: sample λ offspring by bit-flip mutation
      - Selects the best individual among parent ∪ offspring
      - Stops when evaluation budget (budget_factor * n) is exhausted

    Parameters
    ----------
    lam : int
        Number of offspring per iteration (λ ≥ 1).
    budget_factor : int
        Total evaluations budget = budget_factor * n.
    p : Optional[float]
        Mutation rate; defaults to 1/n if None.
    accept_equal : bool
        If True, ties against the parent are accepted (≥ for MAX, ≤ for MIN).
    force_flip_at_least_one : bool
        If True, resamples ell=Bin(n,p) as 1 when ell==0 to ensure at least one bit flips.
    seed : Optional[int]
        RNG seed for reproducibility.

    Tracked attributes (helpful for logging)
    ----------------------------------------
    algorithm_id : int
        Random ID per run (changes on reset()).
    mutation_rate : float
        Effective mutation rate used in the run.
    last_ell : int
        Number of flipped bits in the last sampled offspring (for introspection).
    iterations : int
        Number of iterations (generations) performed so far.
    """
    lam: int = 1
    budget_factor = 5
    seed = 42
    beta: float = 1.5

    _rng = None
    _cdf = None
    def _prepare_powerlaw_sampler(self, n: int):
        # probs ∝ k^{-β} for k=1..n
        k = np.arange(1, n + 1, dtype=float)
        w = k ** (-self.beta)
        p = w / w.sum()
        cdf = np.cumsum(p)
        # store for reuse this run
        self._cdf = cdf

    def _sample_ell(self) -> int:
        # inverse-CDF sampling on cached cdf
        u = self._rng.random()
        return int(1 + np.searchsorted(self._cdf, u))
    
    
    def __call__(self, f, seed=1) -> tuple[float, np.ndarray]:
        """
        Run the (1+λ)EA on a single IOH problem object `f`.

        Returns
        -------
        (fx_best, x_best) : tuple[float, np.ndarray]
            Best fitness found and corresponding solution.
        """
        self.seed = seed
        self._rng = np.random.default_rng(self.seed)

        n = f.meta_data.n_variables
        budget = self.budget_factor * n * n
        self.mutation_rate = 1.0 / n

        self._prepare_powerlaw_sampler(n)
        # Initialize x within bounds (PBO: binary {0,1})
        x = self._rng.integers(0,2,n)
        fx = f(x)


        while f.state.evaluations < budget:
            # Generate λ offspring
            best_y = None
            best_fy = None

            for _ in range(self.lam):
                ell = self._sample_ell()
                while ell < 1:
                    ell = self._sample_ell()
                idx = self._rng.choice(n, size=ell, replace=False)

                y = x.copy()
                # bit flip for {0,1}: new = lb + ub - old (works for 0/1 bounds)
                y[idx] = 1 - y[idx]
                fy = f(y)

                # track last sampled ell (purely for logging/inspection)
                self.last_ell = int(ell)

                if best_y is None:
                    best_y, best_fy = y, fy
                else:
                    if fy >= best_fy:
                        best_y, best_fy = y, fy

            if best_fy >= fx:
                x, fx = best_y, best_fy
            
            if best_fy >= f.optimum.y:
                break

        return fx, x

