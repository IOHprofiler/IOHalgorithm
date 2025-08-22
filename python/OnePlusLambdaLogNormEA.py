import numpy as np
import ioh
import math
# ------------------------ (1+1)EA ------------------------

class OnePlusLambdaLogNormEA:
    lam: int = 10
    budget_factor = 5
    seed = 42
    rng = None

    def sample_conditional_binomial(self, p: float, n: int) -> int:
        """Sample Binomial(n, p), but ensure in [1, n]."""
        k = sum(1 for _ in range(n) if self.rng.random() < p)
        if k < 1:
            k = 1
        elif k > n:
            k = n
        return
    
    def __call__(self, f, seed = 1) -> tuple[float, np.ndarray]:
        """
        Run the (1+λ)EA on a single IOH problem object `f`.

        Returns
        -------
        (fx_best, x_best) : tuple[float, np.ndarray]
            Best fitness found and corresponding solution.
        """
        self.seed = seed
        self.rng = np.random.default_rng(self.seed)

        n = f.meta_data.n_variables
        budget = self.budget_factor * n * n

        # Initialize x within bounds (PBO: binary {0,1})
       
        x = self.rng.integers(0,2,n)
        fx = f(x)

        mutation_rate = 0.2


        while f.state.evaluations < budget:
            # Generate λ offspring
            best_y = None
            best_fy = None
            

            for _ in range(self.lam):
                tmp_mr = 1.0 / (1.0 + (((1.0 - mutation_rate) / mutation_rate) * math.exp(0.22 * rng.normal(0,1))))
                ell = self.sample_conditional_binomial(tmp_mr,n)
                idx = self.rng.choice(n, size=ell, replace=False)

                y = x.copy()
                # bit flip for {0,1}: new = lb + ub - old (works for 0/1 bounds)
                y[idx] = 1 - y[idx]
                fy = f(y)

                if best_y is None:
                    best_y, best_fy, win_r = y, fy, tmp_mr
                else:
                    if fy > best_fy:
                        best_y, best_fy, win_r = y, fy, tmp_mr

    
            if best_fy >= fx:
                x, fx = best_y, best_fy
            
            if best_fy >= f.optimum.y:
                break

            mutation_rate = win_r

        return fx, x

