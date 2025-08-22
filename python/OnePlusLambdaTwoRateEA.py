import numpy as np
import ioh

# ------------------------ (1+1)EA ------------------------

class OnePlusLambdaTwoRateEA:
    lam: int = 10
    budget_factor = 5
    r = 2.0
    seed = 42

    
    def __call__(self, f, seed = 1) -> tuple[float, np.ndarray]:
        """
        Run the (1+λ)EA on a single IOH problem object `f`.

        Returns
        -------
        (fx_best, x_best) : tuple[float, np.ndarray]
            Best fitness found and corresponding solution.
        """
        self.seed = seed
        rng = np.random.default_rng(self.seed)

        n = f.meta_data.n_variables
        budget = self.budget_factor * n * n


        # Initialize x within bounds (PBO: binary {0,1})
       
        x = rng.integers(0,2,n)
        fx = f(x)


        while f.state.evaluations < budget:
            # Generate λ offspring
            best_y = None
            best_fy = None

            p_low = (self.r / (2.0 * n))
            p_high = (2.0 * self.r / n)

            for i in range(self.lam):
                
                if i < self.lam // 2:
                    mutation_rate = p_low
                else:
                    mutation_rate = p_high
                
                ell = rng.binomial(n, mutation_rate)
                while ell < 1:
                    ell = rng.binomial(n, mutation_rate)
                idx = rng.choice(n, size=ell, replace=False)

                y = x.copy()
                # bit flip for {0,1}: new = lb + ub - old (works for 0/1 bounds)
                y[idx] = 1 - y[idx]
                fy = f(y)

                # track last sampled ell (purely for logging/inspection)
                self.last_ell = int(ell)

                if best_y is None:
                    best_y, best_fy, win_mr = y, fy, mutation_rate
                else:
                    if fy >= best_fy:
                        best_y, best_fy, win_mr = y, fy, mutation_rate
                

            if best_fy >= fx:
                x, fx = best_y, best_fy
            
            if best_fy >= f.optimum.y:
                break

            self.r = win_mr * n
            if rng.random() < 0.5:
                if rng.random() < 0.5:
                    self.r = p_low * n
                else:
                    self.r = p_high * n

            # clamp r to [2, n/4]
            if self.r < 2.0:
                self.r = 2.0
            elif self.r > n / 4.0:
                self.r = n / 4.0

        return fx, x
