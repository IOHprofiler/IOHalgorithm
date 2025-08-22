
import numpy as np
import ioh

# ------------------------ (1+1)EA ------------------------

class RLS:
    budget_factor: int = 5       # evaluations = budget_factor * n
    seed = 42


    def reset(self):
        # called between runs by ioh.Experiment; give each run a fresh ID
        self.algorithm_id = np.random.default_rng().integers(10**9)
        self.iterations = 0

    def __call__(self, f, seed= 1):
        """
        Run the algorithm on one IOH problem object (PBO).
        `f` is an IOH function object (callable).
        """
        self.seed = seed
        rng = np.random.default_rng(self.seed)

        n = f.meta_data.n_variables
        budget = self.budget_factor * n * n
        self.mutation_rate = (1.0 / n)
        x = rng.integers(0,2,n)
        fx = f(x)

        while f.state.evaluations < budget:
            # ell = 1
            ell = 1
            idx = rng.choice(n, size=ell, replace=False)

            y = x.copy()
            y[idx] = 1 - y[idx]       # flip bits for 0/1 domain
            fy = f(y)

            
            if fy >= fx:
                x, fx = y, fy

            if fy >= f.optimum.y:
                break
        return fx, x
