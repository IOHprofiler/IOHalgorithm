import numpy as np
import ioh
import math
# ------------------------ (1+1)EA ------------------------

class OnePlusLambdaNormEA:
    lam: int = 10
    budget_factor = 5
    seed = 42
    r : float = 2.0

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

            mean = self.r
            std = math.sqrt(self.r * (1-self.r/n))
            for _ in range(self.lam):
                
                ell = min(math.floor(rng.normal(mean,std)),n / 2)
                while ell < 1:
                    ell = min(math.floor(rng.normal(mean,std)),n / 2)
                idx = rng.choice(n, size=ell, replace=False)

                y = x.copy()
                # bit flip for {0,1}: new = lb + ub - old (works for 0/1 bounds)
                y[idx] = 1 - y[idx]
                fy = f(y)

                if best_y is None:
                    best_y, best_fy, win_r = y, fy, ell
                else:
                    if fy > best_fy:
                        best_y, best_fy, win_r = y, fy, ell

    
            if best_fy >= fx:
                x, fx = best_y, best_fy
            
            if best_fy >= f.optimum.y:
                break

            self.r = win_r

        return fx, x
